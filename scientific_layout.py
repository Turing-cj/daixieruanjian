"""Scientific visual semantics and deterministic layout helpers.

This module deliberately sits between biochemical model construction and the
generic D3 renderer.  It never changes stoichiometry or fluxes; it only
annotates render metadata (roles, positions, routing, labels, and style).
"""
from __future__ import annotations

import math
import itertools
from collections import Counter


GEOMETRY_KEYS = {"x", "y", "label_dx", "label_dy"}


def _lookup(model):
    items = {}
    for kind in ("metabolites", "reactions"):
        for obj in getattr(model, kind):
            items[obj.id] = (kind, obj)
    return items


def _map_info(obj):
    obj.notes = obj.notes or {}
    return obj.notes.setdefault("map_info", {})


def _structure_items(structure):
    """Normalize structure members to inspectable placement specifications."""
    items = []
    for raw in structure.get("items") or []:
        if isinstance(raw, str):
            items.append({"id": raw})
        elif isinstance(raw, dict) and raw.get("id"):
            items.append(dict(raw))
        else:
            raise ValueError(f"invalid scientific_layout structure item: {raw!r}")
    return items


def _place_structures(model, structures):
    """Place ordered pathway skeletons on ellipses / ellipse arcs.

    This layer sits upstream of routing.  It deliberately changes only visual
    geometry and records enough metadata for publication QA to verify that a
    later manual anchor did not silently destroy the intended skeleton.
    """
    lookup = _lookup(model)
    summaries = []
    for index, structure in enumerate(structures or []):
        name = structure.get("name") or f"structure_{index + 1}"
        kind = structure.get("type", "ellipse")
        if kind not in {"ellipse", "ellipse_arc"}:
            raise ValueError(f"unsupported scientific_layout structure type: {kind}")
        center = structure.get("center") or []
        if len(center) != 2:
            raise ValueError(f"{name}: structure center must be [x, y]")
        cx, cy = map(float, center)
        rx = float(structure.get("rx", 0.0))
        ry = float(structure.get("ry", 0.0))
        if rx <= 0 or ry <= 0:
            raise ValueError(f"{name}: structure rx/ry must be positive")
        items = _structure_items(structure)
        if not items:
            continue

        start = float(structure.get("start_angle_deg", 0.0))
        if kind == "ellipse_arc":
            if "end_angle_deg" not in structure:
                raise ValueError(f"{name}: ellipse_arc requires end_angle_deg")
            end = float(structure["end_angle_deg"])
        else:
            direction = str(structure.get("direction", "clockwise")).lower()
            sweep = 360.0 if direction != "counterclockwise" else -360.0
            # Do not duplicate the first point at 360 degrees for a closed ring.
            end = start + sweep * (len(items) - 1) / max(len(items), 1)

        if len(items) == 1:
            generated_angles = [start]
        else:
            generated_angles = [
                start + (end - start) * item_index / (len(items) - 1)
                for item_index in range(len(items))
            ]

        placed = []
        for order_index, (item, generated_angle) in enumerate(zip(items, generated_angles)):
            item_id = item["id"]
            if item_id not in lookup:
                raise ValueError(f"{name}: structure references unknown id: {item_id}")
            angle = float(item.get("angle_deg", generated_angle))
            theta = math.radians(angle)
            base_dx, base_dy = rx * math.cos(theta), ry * math.sin(theta)
            radial_length = math.hypot(base_dx, base_dy) or 1.0
            radial_offset = float(item.get("radial_offset", 0.0))
            x = cx + base_dx + base_dx / radial_length * radial_offset
            y = cy + base_dy + base_dy / radial_length * radial_offset
            _, obj = lookup[item_id]
            info = _map_info(obj)
            info["x"], info["y"] = x, y
            info["layout_structure"] = name
            info["layout_structure_type"] = kind
            info["layout_structure_center"] = [cx, cy]
            info["layout_structure_radius"] = [rx, ry]
            info["layout_theta_deg"] = angle
            info["layout_radial_offset"] = radial_offset
            info["layout_structure_order_index"] = order_index
            if structure.get("region"):
                info["layout_region"] = structure["region"]
            placed.append(item_id)

        summaries.append({
            "name": name,
            "type": kind,
            "items": placed,
            "center": [cx, cy],
            "rx": rx,
            "ry": ry,
            "start_angle_deg": start,
            "end_angle_deg": end,
        })
    return summaries


def _optimize_structure_orders(model, structures, route_policy):
    """Exact small-N angular-order search for circular pathway structures.

    Node order is a higher-level variable than spline routing.  Searching it
    before local route repair removes topological crossings instead of asking
    waypoints to route around a bad embedding.  Only structures that opt in
    are touched, and the search is capped so the behavior stays transparent.
    """
    lookup = _lookup(model)
    channel_spacing = float(route_policy.get("channel_spacing", 22.0))
    obstacle_clearance = float(route_policy.get("obstacle_clearance", 42.0))
    default_limit = int(route_policy.get("max_exact_structure_order_size", 7))
    results = []

    for structure in structures or []:
        if not structure.get("optimize_order"):
            continue
        items = _structure_items(structure)
        name = structure.get("name", "unnamed")
        limit = int(structure.get("max_exact_order_size", default_limit))
        if len(items) <= 1 or len(items) > limit:
            results.append({
                "name": name,
                "optimized": False,
                "reason": "structure size outside exact-search range",
            })
            continue

        item_ids = [item["id"] for item in items]
        if any(item_id not in lookup for item_id in item_ids):
            raise ValueError(f"{name}: optimize_order references unknown item")

        # Slots are geometry, not identities.  Keep the original angular/radial
        # slots and search which biochemical node occupies each slot.
        slots = []
        for item_id in item_ids:
            info = _map_info(lookup[item_id][1])
            slots.append({
                "x": float(info["x"]),
                "y": float(info["y"]),
                "theta": float(info["layout_theta_deg"]),
                "radial_offset": float(info.get("layout_radial_offset", 0.0)),
            })

        before_edges = _all_compiled_corridors(model, channel_spacing)
        before_quality = _route_defect_score(model, before_edges, obstacle_clearance)
        before_length = _route_length(before_edges)
        original_order = tuple(item_ids)
        best = None

        for permutation in itertools.permutations(item_ids):
            for slot_index, item_id in enumerate(permutation):
                info = _map_info(lookup[item_id][1])
                slot = slots[slot_index]
                info["x"], info["y"] = slot["x"], slot["y"]
                info["layout_theta_deg"] = slot["theta"]
                info["layout_radial_offset"] = slot["radial_offset"]
                info["layout_structure_order_index"] = slot_index

            candidate_edges = _all_compiled_corridors(model, channel_spacing)
            quality = _route_defect_score(model, candidate_edges, obstacle_clearance)
            route_length = _route_length(candidate_edges)
            displacement = sum(
                int(item_id != original_order[index])
                for index, item_id in enumerate(permutation)
            )
            score = quality + (route_length, displacement, permutation)
            if best is None or score < best[0]:
                best = (score, permutation)

        _, best_order = best
        for slot_index, item_id in enumerate(best_order):
            info = _map_info(lookup[item_id][1])
            slot = slots[slot_index]
            info["x"], info["y"] = slot["x"], slot["y"]
            info["layout_theta_deg"] = slot["theta"]
            info["layout_radial_offset"] = slot["radial_offset"]
            info["layout_structure_order_index"] = slot_index

        after_edges = _all_compiled_corridors(model, channel_spacing)
        after_quality = _route_defect_score(model, after_edges, obstacle_clearance)
        results.append({
            "name": name,
            "optimized": True,
            "order_before": list(original_order),
            "order_after": list(best_order),
            "defects_before": before_quality[0],
            "defects_after": after_quality[0],
            "crossings_before": before_quality[1],
            "crossings_after": after_quality[1],
            "obstacles_before": before_quality[2],
            "obstacles_after": after_quality[2],
            "route_length_before": before_length,
            "route_length_after": _route_length(after_edges),
        })
    return results


def apply_scientific_layout(model, layout):
    """Compile optional scientific-layout semantics into D3Flux map_info.

    Supported primitives intentionally stay small and inspectable:

    * roles: semantic role -> list of metabolite/reaction IDs
    * lanes: evenly place ordered IDs along an x or y lane
    * anchors: explicit final positions for exceptional junctions/cycles
    * metabolite_styles / reaction_styles: visual metadata overrides
    * role_palette: optional semantic-role colors

    The layout file remains the source of truth; no biochemical state changes.
    """
    if layout.get('circular_backbone', {}).get('enabled'):
        from circular_layout import compile_circular_layout
        return compile_circular_layout(model, layout)
    spec = layout.get("scientific_layout") or {}
    if not spec.get("enabled"):
        return {"enabled": False, "role_counts": {}, "lanes": 0, "anchors": 0}

    lookup = _lookup(model)
    role_palette = spec.get("role_palette", {})
    role_counts = Counter()

    for role, ids in spec.get("roles", {}).items():
        for item_id in ids:
            if item_id not in lookup:
                raise ValueError(f"scientific_layout role references unknown id: {item_id}")
            _, obj = lookup[item_id]
            info = _map_info(obj)
            info["scientific_role"] = role
            if role in role_palette:
                info["color"] = role_palette[role]
            role_counts[role] += 1

    for lane in spec.get("lanes", []):
        ids = lane.get("items") or []
        if not ids:
            continue
        axis = lane.get("axis", "x")
        if axis not in {"x", "y"}:
            raise ValueError(f"scientific_layout lane axis must be x/y, got: {axis}")
        start = float(lane["start"])
        end = float(lane["end"])
        fixed = float(lane["fixed"])
        step = 0.0 if len(ids) == 1 else (end - start) / (len(ids) - 1)
        for index, item_id in enumerate(ids):
            if item_id not in lookup:
                raise ValueError(f"scientific_layout lane references unknown id: {item_id}")
            _, obj = lookup[item_id]
            info = _map_info(obj)
            variable = start + step * index
            if axis == "x":
                info["x"], info["y"] = variable, fixed
            else:
                info["x"], info["y"] = fixed, variable
            info["layout_lane"] = lane.get("name", "unnamed")

    structure_summary = _place_structures(model, spec.get("structures", []))

    for item_id, position in spec.get("anchors", {}).items():
        if item_id not in lookup:
            raise ValueError(f"scientific_layout anchor references unknown id: {item_id}")
        _, obj = lookup[item_id]
        info = _map_info(obj)
        info["x"] = float(position["x"])
        info["y"] = float(position["y"])

    for key, styles in (("metabolites", spec.get("metabolite_styles", {})),
                        ("reactions", spec.get("reaction_styles", {}))):
        expected_kind = key
        for item_id, values in styles.items():
            if item_id not in lookup:
                raise ValueError(f"scientific_layout style references unknown id: {item_id}")
            kind, obj = lookup[item_id]
            if kind != expected_kind:
                raise ValueError(f"scientific_layout {key} style applied to {kind[:-1]}: {item_id}")
            _map_info(obj).update(values)

    route_policy = spec.get("route_policy") or {}
    structure_order_summary = _optimize_structure_orders(
        model, spec.get("structures", []), route_policy
    )
    route_summary = _compile_route_grammar(model, route_policy)

    label_policy = spec.get("label_policy") or {}
    label_summary = _resolve_label_offsets(model, label_policy)
    summary = {
        "enabled": True,
        "role_counts": dict(sorted(role_counts.items())),
        "lanes": len(spec.get("lanes", [])),
        "structures": structure_summary,
        "structure_order_optimization": structure_order_summary,
        "anchors": len(spec.get("anchors", {})),
        "routes": route_summary,
        "labels": label_summary,
    }
    model.notes.setdefault("scientific_layout", {}).update(summary)
    return summary


def _unit(dx, dy):
    length = math.hypot(dx, dy) or 1.0
    return dx / length, dy / length


def _cardinal_port(dx, dy):
    """Choose a stable N/S/E/W node port that faces the next waypoint."""
    if abs(dx) >= abs(dy):
        return "E" if dx >= 0 else "W"
    return "S" if dy >= 0 else "N"


def _visible_reaction_sides(reaction):
    """Return visible non-local reactants/products for route semantics.

    A biochemical reaction is a hyperedge. The legacy D3Flux renderer later
    pairs reactants/products only as a drawing convenience, so routing code
    must keep the two reaction sides available independently instead of
    treating those arbitrary pairs as biochemical edges.
    """
    info = _map_info(reaction)
    local = set((info.get("cofactors") or {}).keys())
    reactants = [m for m in reaction.reactants if m.id not in local and not _map_info(m).get("hidden")]
    products = [m for m in reaction.products if m.id not in local and not _map_info(m).get("hidden")]
    return reactants, products


def _edge_pairs(reaction):
    """Mirror D3Flux's deterministic drawing pairs.

    The pairs are not biochemical substrate-to-product assignments. They only
    ensure every visible endpoint participates in at least one rendered path.
    Multi-endpoint reactions can therefore use shared reaction-bus geometry
    while retaining D3Flux's existing bilink data model.
    """
    reactants, products = _visible_reaction_sides(reaction)
    if not reactants or not products:
        return []
    if len(reactants) >= len(products):
        return [(reactant, products[index % len(products)])
                for index, reactant in enumerate(reactants)]
    return [(reactants[index % len(reactants)], product)
            for index, product in enumerate(products)]


def _reaction_bus_geometry(reaction, half_length=34.0):
    """Build a compact hyperedge spine for a multi-endpoint reaction.

    All reactant branches converge on the input join and all product branches
    leave from the output join. The short join-to-join segment is the shared
    reaction trunk. Its axis follows reactant-centroid to product-centroid, so
    arbitrary substrate/product drawing pairs never become visible as long
    crossing curves.
    """
    info = _map_info(reaction)
    reactants, products = _visible_reaction_sides(reaction)
    if not reactants or not products or not {"x", "y"} <= info.keys():
        return None

    def centroid(items):
        points = [
            (float(_map_info(item)["x"]), float(_map_info(item)["y"]))
            for item in items
            if {"x", "y"} <= _map_info(item).keys()
        ]
        if not points:
            return None
        return (
            sum(point[0] for point in points) / len(points),
            sum(point[1] for point in points) / len(points),
        )

    reactant_center = centroid(reactants)
    product_center = centroid(products)
    if reactant_center is None or product_center is None:
        return None

    same_centroid = (
        math.isclose(product_center[0], reactant_center[0], abs_tol=1e-9)
        and math.isclose(product_center[1], reactant_center[1], abs_tol=1e-9)
    )
    if same_centroid:
        rx, ry = float(info["x"]), float(info["y"])
        same_as_hub = (
            math.isclose(product_center[0], rx, abs_tol=1e-9)
            and math.isclose(product_center[1], ry, abs_tol=1e-9)
        )
        if same_as_hub:
            ux, uy = 1.0, 0.0
        else:
            ux, uy = _unit(product_center[0] - rx, product_center[1] - ry)
    else:
        ux, uy = _unit(
            product_center[0] - reactant_center[0],
            product_center[1] - reactant_center[1],
        )

    half = max(12.0, float(half_length))
    rx, ry = float(info["x"]), float(info["y"])
    return {
        "axis": [ux, uy],
        "input": [rx - ux * half, ry - uy * half],
        "output": [rx + ux * half, ry + uy * half],
        "half_length": half,
        "reactants": [item.id for item in reactants],
        "products": [item.id for item in products],
    }


def _segments_intersect(a, b, c, d):
    """Proper segment intersection; touching endpoints are not crossings."""
    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    if a == c or a == d or b == c or b == d:
        return False
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    return (o1 * o2 < 0) and (o3 * o4 < 0)


def _point_segment_distance(point, a, b):
    vx, vy = b[0] - a[0], b[1] - a[1]
    wx, wy = point[0] - a[0], point[1] - a[1]
    length2 = vx * vx + vy * vy
    if length2 <= 1e-12:
        return math.hypot(wx, wy)
    t = max(0.0, min(1.0, (wx * vx + wy * vy) / length2))
    px, py = a[0] + t * vx, a[1] + t * vy
    return math.hypot(point[0] - px, point[1] - py)


def _segment_distance(a, b, c, d):
    if _segments_intersect(a, b, c, d):
        return 0.0
    return min(
        _point_segment_distance(a, c, d),
        _point_segment_distance(b, c, d),
        _point_segment_distance(c, a, b),
        _point_segment_distance(d, a, b),
    )


def _polyline_distance(first, second):
    best = float("inf")
    for a, b in zip(first["points"], first["points"][1:]):
        for c, d in zip(second["points"], second["points"][1:]):
            best = min(best, _segment_distance(a, b, c, d))
    return best


def _route_visual_metrics(model, channel_spacing, edge_clearance=18.0,
                          max_route_stretch=2.4):
    """Audit visual clutter that strict segment crossings do not capture.

    The metric is intentionally route-level and conservative. Edges that
    belong to the same reaction or share a biochemical endpoint are excluded
    from near-edge penalties because their local convergence is intentional.
    """
    compiled = _all_compiled_corridors(model, channel_spacing)
    near_pairs = []
    for index, first in enumerate(compiled):
        for second in compiled[index + 1:]:
            if first["reaction"] == second["reaction"]:
                continue
            if {first["source"], first["target"]} & {second["source"], second["target"]}:
                continue
            distance = _polyline_distance(first, second)
            if 1e-6 < distance < edge_clearance:
                near_pairs.append({
                    "first": f'{first["reaction"]}:{first["source"]}->{first["target"]}',
                    "second": f'{second["reaction"]}:{second["source"]}->{second["target"]}',
                    "distance": distance,
                })

    stretch = []
    reaction_positions = {
        reaction.id: (float(_map_info(reaction)["x"]), float(_map_info(reaction)["y"]))
        for reaction in model.reactions
        if {"x", "y"} <= _map_info(reaction).keys()
    }
    for edge in compiled:
        direct = math.hypot(
            edge["points"][-1][0] - edge["points"][0][0],
            edge["points"][-1][1] - edge["points"][0][1],
        )
        reaction = None
        try:
            reaction = model.reactions.get_by_id(edge["reaction"])
        except KeyError:
            pass
        if reaction is not None and _map_info(reaction).get("route_style") == "reaction_bus":
            # The hyperedge bus is itself the semantic baseline.  Do not
            # compare it against a source->single-hub->target shortcut, which
            # would reward collapsing a grouped reaction back into bilinks.
            semantic_baseline = sum(
                math.hypot(b[0] - a[0], b[1] - a[1])
                for a, b in zip(edge["points"], edge["points"][1:])
            )
            hub = None
        else:
            hub = reaction_positions.get(edge["reaction"])
        if hub:
            semantic_baseline = (
                math.hypot(hub[0] - edge["points"][0][0],
                           hub[1] - edge["points"][0][1])
                + math.hypot(edge["points"][-1][0] - hub[0],
                             edge["points"][-1][1] - hub[1])
            )
        elif reaction is None or _map_info(reaction).get("route_style") != "reaction_bus":
            semantic_baseline = direct
        routed = sum(
            math.hypot(b[0] - a[0], b[1] - a[1])
            for a, b in zip(edge["points"], edge["points"][1:])
        )
        ratio = routed / max(semantic_baseline, 1e-9)
        if ratio > max_route_stretch:
            stretch.append({
                "edge": f'{edge["reaction"]}:{edge["source"]}->{edge["target"]}',
                "ratio": ratio,
                "routed_length": routed,
                "direct_length": direct,
                "semantic_baseline_length": semantic_baseline,
            })

    return {
        "edge_clearance_threshold": edge_clearance,
        "near_edge_pairs": near_pairs,
        "near_edge_count": len(near_pairs),
        "max_route_stretch": max_route_stretch,
        "route_stretch_violations": stretch,
        "route_stretch_violation_count": len(stretch),
    }


def _edge_obstacle_hits(model, edge, clearance=42.0):
    """Return unrelated visible nodes approached too closely by one route."""
    excluded = {edge["reaction"], edge["source"], edge["target"]}
    hits = []
    for kind in ("metabolites", "reactions"):
        for obj in getattr(model, kind):
            if obj.id in excluded:
                continue
            info = _map_info(obj)
            if info.get("hidden") or not {"x", "y"} <= info.keys():
                continue
            point = (float(info["x"]), float(info["y"]))
            if any(_point_segment_distance(point, a, b) < clearance
                   for a, b in zip(edge["points"], edge["points"][1:])):
                hits.append(obj.id)
    return sorted(set(hits))


def _route_obstacle_hits(model, compiled_edges, clearance=42.0):
    hits = []
    for edge in compiled_edges:
        obstacles = _edge_obstacle_hits(model, edge, clearance)
        if obstacles:
            hits.append({
                "edge": f'{edge["reaction"]}:{edge["source"]}->{edge["target"]}',
                "obstacles": obstacles,
            })
    return hits


def _route_defect_score(model, compiled_edges, clearance=42.0):
    crossings = _estimate_route_crossings(compiled_edges)
    obstacle_items = _route_obstacle_hits(model, compiled_edges, clearance)
    obstacle_count = sum(len(item["obstacles"]) for item in obstacle_items)
    return crossings + obstacle_count, crossings, obstacle_count


def _route_length(compiled_edges):
    return sum(
        sum(math.hypot(b[0] - a[0], b[1] - a[1])
            for a, b in zip(edge["points"], edge["points"][1:]))
        for edge in compiled_edges
    )


def _route_compaction_objective(model, compiled_edges, channel_spacing,
                                flexible_origins):
    """Prefer shorter zero-defect routes without large visual displacements."""
    node_move = 0.0
    for metabolite in model.metabolites:
        if metabolite.id not in flexible_origins:
            continue
        info = _map_info(metabolite)
        ox, oy = flexible_origins[metabolite.id]
        node_move += math.hypot(float(info["x"]) - ox, float(info["y"]) - oy)

    offset_move = 0.0
    waypoint_detour = 0.0
    for reaction in model.reactions:
        info = _map_info(reaction)
        spacing = float(info.get("route_channel_spacing", channel_spacing))
        reaction_channel = int(info.get("route_channel", 0))
        overrides = info.get("edge_route_overrides") or {}
        pairs = {f"{source.id}->{target.id}": (source, target)
                 for source, target in _edge_pairs(reaction)}
        for key, override in overrides.items():
            pair = pairs.get(key)
            if not pair:
                continue
            source, target = pair
            style = override.get("style", info.get("route_style", "bezier"))
            if style == "arc":
                style = "bundle_spline"
            baseline_channel = int(override.get("channel", reaction_channel))
            baseline_offset = baseline_channel * spacing if style == "bundle_spline" else 0.0
            edge_offset = float(override.get("channel_offset", baseline_offset))
            if "channel_offset" in override:
                offset_move += abs(edge_offset - baseline_offset)

            sx, sy = float(_map_info(source)["x"]), float(_map_info(source)["y"])
            tx, ty = float(_map_info(target)["x"]), float(_map_info(target)["y"])
            rx, ry = float(info["x"]), float(info["y"])
            ux, uy = _unit(tx - sx, ty - sy)
            hx, hy = rx - uy * edge_offset, ry + ux * edge_offset
            for field, a, b in (
                ("waypoint_before", (sx, sy), (hx, hy)),
                ("waypoint_after", (hx, hy), (tx, ty)),
            ):
                waypoint = override.get(field)
                if not waypoint:
                    continue
                mx, my = (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0
                waypoint_detour += math.hypot(
                    float(waypoint["x"]) - mx, float(waypoint["y"]) - my
                )

    length = _route_length(compiled_edges)
    score = length + 0.35 * node_move + 0.15 * waypoint_detour + 0.10 * offset_move
    return score, length, node_move, waypoint_detour, offset_move


def _compact_zero_defect_routes(model, channel_spacing, flexible_origins,
                                obstacle_clearance=42.0, max_node_nudge=160,
                                node_step=40, max_edge_channel=6,
                                max_passes=4):
    """Shorten an already valid route while preserving zero route defects.

    This is deliberately a second phase. Correctness is frozen first
    (crossings/obstacles == 0); only then may geometry be shortened.
    """
    compiled = _all_compiled_corridors(model, channel_spacing)
    defects = _route_defect_score(model, compiled, obstacle_clearance)
    before = _route_compaction_objective(
        model, compiled, channel_spacing, flexible_origins
    )
    if defects[0] != 0:
        return {
            "enabled": True,
            "applied": False,
            "reason": "route still has defects",
            "objective_before": before,
            "objective_after": before,
            "changes": [],
        }

    changes = []

    def current_objective():
        edges = _all_compiled_corridors(model, channel_spacing)
        if _route_defect_score(model, edges, obstacle_clearance)[0] != 0:
            return None
        return _route_compaction_objective(
            model, edges, channel_spacing, flexible_origins
        )

    for _ in range(max_passes):
        changed = False

        # 1) Flexible nodes: allow a slightly larger visual move if the total
        # route becomes materially shorter while defects remain zero.
        for metabolite in sorted(
                (m for m in model.metabolites if m.id in flexible_origins),
                key=lambda item: item.id):
            info = _map_info(metabolite)
            current = (float(info["x"]), float(info["y"]))
            ox, oy = flexible_origins[metabolite.id]
            baseline = current_objective()
            best = (baseline, current)
            for dx in range(-max_node_nudge, max_node_nudge + 1, node_step):
                for dy in range(-max_node_nudge, max_node_nudge + 1, node_step):
                    candidate = (ox + dx, oy + dy)
                    info["x"], info["y"] = candidate
                    objective = current_objective()
                    if objective is not None and objective < best[0]:
                        best = (objective, candidate)
            info["x"], info["y"] = best[1]
            if best[1] != current:
                changes.append({
                    "type": "node",
                    "id": metabolite.id,
                    "from": list(current),
                    "to": list(best[1]),
                })
                changed = True

        # 2) Per-edge channel offsets: pull them back toward their semantic
        # bundle channel when the zero-defect invariant permits it.
        for reaction in model.reactions:
            info = _map_info(reaction)
            overrides = info.get("edge_route_overrides") or {}
            spacing = float(info.get("route_channel_spacing", channel_spacing))
            base_channel = int(info.get("route_channel", 0))
            for key, override in sorted(overrides.items()):
                if "channel_offset" not in override:
                    continue
                current = float(override["channel_offset"])
                baseline = current_objective()
                best = (baseline, current)
                values = [float(index * spacing)
                          for index in range(-max_edge_channel, max_edge_channel + 1)]
                for value in values:
                    override["channel_offset"] = value
                    objective = current_objective()
                    if objective is not None and objective < best[0]:
                        best = (objective, value)
                override["channel_offset"] = best[1]
                if not math.isclose(best[1], current, abs_tol=1e-9):
                    changes.append({
                        "type": "channel_offset",
                        "reaction": reaction.id,
                        "edge": key,
                        "from": current,
                        "to": best[1],
                    })
                    changed = True

        # 3) Waypoints: contract each detour toward the midpoint of the
        # underlying segment. Shared-knot cubic rendering keeps the result
        # smooth as the detour contracts.
        for reaction in model.reactions:
            info = _map_info(reaction)
            overrides = info.get("edge_route_overrides") or {}
            pairs = {f"{source.id}->{target.id}": (source, target)
                     for source, target in _edge_pairs(reaction)}
            for key, override in sorted(overrides.items()):
                pair = pairs.get(key)
                if not pair:
                    continue
                source, target = pair
                spacing = float(info.get("route_channel_spacing", channel_spacing))
                style = override.get("style", info.get("route_style", "bezier"))
                if style == "arc":
                    style = "bundle_spline"
                channel = int(override.get("channel", info.get("route_channel", 0)))
                default_offset = channel * spacing if style == "bundle_spline" else 0.0
                edge_offset = float(override.get("channel_offset", default_offset))
                sx, sy = float(_map_info(source)["x"]), float(_map_info(source)["y"])
                tx, ty = float(_map_info(target)["x"]), float(_map_info(target)["y"])
                rx, ry = float(info["x"]), float(info["y"])
                ux, uy = _unit(tx - sx, ty - sy)
                hx, hy = rx - uy * edge_offset, ry + ux * edge_offset
                for field, a, b in (
                    ("waypoint_before", (sx, sy), (hx, hy)),
                    ("waypoint_after", (hx, hy), (tx, ty)),
                ):
                    if field not in override:
                        continue
                    current = dict(override[field])
                    mx, my = (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0
                    baseline = current_objective()
                    best = (baseline, current)
                    for alpha in (0.0, 0.2, 0.3, 0.4, 0.5, 0.6,
                                  0.7, 0.8, 0.9, 1.0):
                        candidate = {
                            "x": mx + (current["x"] - mx) * alpha,
                            "y": my + (current["y"] - my) * alpha,
                        }
                        override[field] = candidate
                        objective = current_objective()
                        if objective is not None and objective < best[0]:
                            best = (objective, dict(candidate))
                    override[field] = best[1]
                    if (not math.isclose(best[1]["x"], current["x"], abs_tol=1e-9)
                            or not math.isclose(best[1]["y"], current["y"], abs_tol=1e-9)):
                        changes.append({
                            "type": "waypoint",
                            "reaction": reaction.id,
                            "edge": key,
                            "field": field,
                            "from": current,
                            "to": dict(best[1]),
                        })
                        changed = True

        if not changed:
            break

    final_edges = _all_compiled_corridors(model, channel_spacing)
    after = _route_compaction_objective(
        model, final_edges, channel_spacing, flexible_origins
    )
    return {
        "enabled": True,
        "applied": bool(changes),
        "objective_before": before,
        "objective_after": after,
        "route_length_before": before[1],
        "route_length_after": after[1],
        "changes": changes,
    }


def _raw_route_crossings(compiled_edges):
    """Return intersecting edge records, excluding intentional shared nodes."""
    crossings = []
    for index, first in enumerate(compiled_edges):
        for second in compiled_edges[index + 1:]:
            if first["reaction"] == second["reaction"]:
                continue
            if {first["source"], first["target"]} & {second["source"], second["target"]}:
                continue
            found = False
            for first_index, (a, b) in enumerate(zip(first["points"], first["points"][1:])):
                for second_index, (c, d) in enumerate(zip(second["points"], second["points"][1:])):
                    if _segments_intersect(a, b, c, d):
                        crossings.append((first, second, first_index, second_index))
                        found = True
                        break
                if found:
                    break
    return crossings


def _route_crossing_pairs(compiled_edges):
    """Return serializable inter-reaction corridor crossings for route QA."""
    return [{
        "first": f'{first["reaction"]}:{first["source"]}->{first["target"]}',
        "second": f'{second["reaction"]}:{second["source"]}->{second["target"]}',
    } for first, second, _, _ in _raw_route_crossings(compiled_edges)]


def _estimate_route_crossings(compiled_edges):
    """Cheap route audit using source→shifted-hub→target polylines."""
    return len(_route_crossing_pairs(compiled_edges))


def _all_compiled_corridors(model, channel_spacing):
    compiled = []
    for reaction in model.reactions:
        compiled.extend(_compiled_corridors_for_reaction(reaction, channel_spacing))
    return compiled


def _effective_edge_offset(reaction, edge_key, channel_spacing):
    info = _map_info(reaction)
    override = (info.get("edge_route_overrides") or {}).get(edge_key, {})
    style = override.get("style", info.get("route_style", "bezier"))
    if style == "arc":
        style = "bundle_spline"
    if "channel_offset" in override:
        return float(override["channel_offset"])
    spacing = float(info.get("route_channel_spacing", channel_spacing))
    channel = int(override.get("channel", info.get("route_channel", 0)))
    return channel * spacing if style == "bundle_spline" else 0.0


def _optimize_edge_routes(model, channel_spacing, max_channel=6, max_iterations=12,
                          obstacle_clearance=42.0):
    """Resolve residual crossings by moving only the conflicting edge corridor.

    Reaction-level channels cannot solve every topology because one reaction
    can contain multiple source→target bilinks.  This second stage therefore
    searches per-edge channel offsets, preserving the reaction hub, ports,
    stoichiometry, and all non-conflicting edges.
    """
    by_id = {reaction.id: reaction for reaction in model.reactions}
    before_edges = _all_compiled_corridors(model, channel_spacing)
    before_quality = _route_defect_score(model, before_edges, obstacle_clearance)
    changes = []

    for _ in range(max_iterations):
        compiled = _all_compiled_corridors(model, channel_spacing)
        raw = _raw_route_crossings(compiled)
        if not raw:
            break
        current_quality = _route_defect_score(model, compiled, obstacle_clearance)
        candidate_edges = {}
        for first, second, _, _ in raw:
            for edge in (first, second):
                reaction = by_id[edge["reaction"]]
                if _map_info(reaction).get("route_style") == "reaction_bus":
                    # Hyperedge branches share one reaction spine. Moving one
                    # arbitrary D3Flux bilink independently would recreate the
                    # spider-web failure this representation is meant to avoid.
                    continue
                key = (edge["reaction"], f'{edge["source"]}->{edge["target"]}')
                candidate_edges[key] = edge

        best = None
        offsets = [0.0]
        for index in range(1, max_channel + 1):
            offsets.extend([-index * channel_spacing, index * channel_spacing])

        for (reaction_id, edge_key), _edge in sorted(candidate_edges.items()):
            reaction = by_id[reaction_id]
            info = _map_info(reaction)
            overrides = info.setdefault("edge_route_overrides", {})
            override = overrides.setdefault(edge_key, {})
            if override.get("locked"):
                continue
            had_offset = "channel_offset" in override
            previous_offset = override.get("channel_offset")
            effective_before = _effective_edge_offset(reaction, edge_key, channel_spacing)

            for offset in offsets:
                if math.isclose(offset, effective_before, abs_tol=1e-9):
                    continue
                override["channel_offset"] = float(offset)
                candidate = _all_compiled_corridors(model, channel_spacing)
                quality = _route_defect_score(model, candidate, obstacle_clearance)
                movement = abs(offset - effective_before)
                objective = quality + (movement, abs(offset), reaction_id, edge_key, offset)
                if best is None or objective < best[0]:
                    best = (objective, reaction, edge_key, float(offset),
                            effective_before, had_offset, previous_offset)

            if had_offset:
                override["channel_offset"] = previous_offset
            else:
                override.pop("channel_offset", None)

        if best is None or best[0][:3] >= current_quality:
            break

        _, reaction, edge_key, offset, effective_before, _, _ = best
        override = _map_info(reaction).setdefault(
            "edge_route_overrides", {}
        ).setdefault(edge_key, {})
        override["channel_offset"] = offset
        changes.append({
            "reaction": reaction.id,
            "edge": edge_key,
            "offset_before": effective_before,
            "offset_after": offset,
        })

    final_edges = _all_compiled_corridors(model, channel_spacing)
    final_quality = _route_defect_score(model, final_edges, obstacle_clearance)
    return {
        "enabled": True,
        "defects_before": before_quality[0],
        "defects_after": final_quality[0],
        "crossings_before": before_quality[1],
        "crossings_after": final_quality[1],
        "obstacles_before": before_quality[2],
        "obstacles_after": final_quality[2],
        "changes": changes,
    }


def _optimize_route_hubs(model, channel_spacing, max_nudge=160, step=40,
                         max_iterations=6, obstacle_clearance=42.0):
    """Last-resort local topology repair for explicitly flexible reaction hubs.

    Metabolite coordinates stay fixed. Only reactions carrying
    route_hub_flexible may move, and a move is committed only when it
    strictly reduces the measured route defect tuple.
    """
    by_id = {reaction.id: reaction for reaction in model.reactions}
    before_edges = _all_compiled_corridors(model, channel_spacing)
    before_quality = _route_defect_score(model, before_edges, obstacle_clearance)
    changes = []

    offsets = [(dx, dy)
               for dx in range(-max_nudge, max_nudge + 1, step)
               for dy in range(-max_nudge, max_nudge + 1, step)
               if dx or dy]
    offsets.sort(key=lambda value: (math.hypot(*value), abs(value[0]) + abs(value[1]),
                                    value[0], value[1]))

    for _ in range(max_iterations):
        compiled = _all_compiled_corridors(model, channel_spacing)
        raw = _raw_route_crossings(compiled)
        current_quality = _route_defect_score(model, compiled, obstacle_clearance)
        obstacle_hits = _route_obstacle_hits(model, compiled, obstacle_clearance)
        if current_quality[0] == 0:
            break
        candidates = set()
        for first, second, _, _ in raw:
            for edge in (first, second):
                reaction = by_id[edge["reaction"]]
                if _map_info(reaction).get("route_hub_flexible"):
                    candidates.add(reaction.id)
        for item in obstacle_hits:
            reaction_id = item["edge"].split(":", 1)[0]
            reaction = by_id.get(reaction_id)
            if reaction is not None and _map_info(reaction).get("route_hub_flexible"):
                candidates.add(reaction_id)
        if not candidates:
            break

        best = None
        for reaction_id in sorted(candidates):
            reaction = by_id[reaction_id]
            info = _map_info(reaction)
            origin_x, origin_y = float(info["x"]), float(info["y"])
            for dx, dy in offsets:
                info["x"], info["y"] = origin_x + dx, origin_y + dy
                candidate = _all_compiled_corridors(model, channel_spacing)
                quality = _route_defect_score(model, candidate, obstacle_clearance)
                objective = quality + (math.hypot(dx, dy), abs(dx) + abs(dy),
                             reaction_id, dx, dy)
                if best is None or objective < best[0]:
                    best = (objective, reaction, origin_x, origin_y, dx, dy)
            info["x"], info["y"] = origin_x, origin_y

        if best is None or best[0][:3] >= current_quality:
            break
        _, reaction, origin_x, origin_y, dx, dy = best
        info = _map_info(reaction)
        info["x"], info["y"] = origin_x + dx, origin_y + dy
        changes.append({
            "reaction": reaction.id,
            "from": [origin_x, origin_y],
            "to": [info["x"], info["y"]],
        })

    final = _all_compiled_corridors(model, channel_spacing)
    final_quality = _route_defect_score(model, final, obstacle_clearance)
    return {
        "enabled": True,
        "defects_before": before_quality[0],
        "defects_after": final_quality[0],
        "crossings_before": before_quality[1],
        "crossings_after": final_quality[1],
        "obstacles_before": before_quality[2],
        "obstacles_after": final_quality[2],
        "changes": changes,
    }


def _optimize_route_waypoints(model, channel_spacing, max_detour=240,
                              detour_step=40, max_iterations=8,
                              obstacle_clearance=42.0):
    """Insert local detour waypoints for residual crossing segments.

    The search only touches the specific segment that participates in a
    crossing. Metabolite and reaction-hub coordinates remain unchanged.
    """
    by_id = {reaction.id: reaction for reaction in model.reactions}
    before_edges = _all_compiled_corridors(model, channel_spacing)
    before_quality = _route_defect_score(model, before_edges, obstacle_clearance)
    changes = []
    distances = []
    for distance in range(detour_step, max_detour + 1, detour_step):
        distances.extend([-float(distance), float(distance)])

    for _ in range(max_iterations):
        compiled = _all_compiled_corridors(model, channel_spacing)
        raw = _raw_route_crossings(compiled)
        if not raw:
            break
        current_quality = _route_defect_score(model, compiled, obstacle_clearance)
        best = None

        for first, second, first_segment, second_segment in raw:
            for edge, segment_index in ((first, first_segment), (second, second_segment)):
                reaction = by_id[edge["reaction"]]
                info = _map_info(reaction)
                if info.get("route_style") == "reaction_bus":
                    # A bus is a grouped hyperedge; detouring one paired
                    # bilink would split the shared trunk into pair-specific
                    # long curves.
                    continue
                edge_key = f'{edge["source"]}->{edge["target"]}'
                override = info.setdefault("edge_route_overrides", {}).setdefault(edge_key, {})
                if override.get("locked"):
                    continue

                hub_index = 1 + int(bool(override.get("waypoint_before")))
                field = "waypoint_before" if segment_index < hub_index else "waypoint_after"
                if field in override:
                    continue

                a = edge["points"][segment_index]
                b = edge["points"][segment_index + 1]
                mx, my = (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0
                ux, uy = _unit(b[0] - a[0], b[1] - a[1])
                nx, ny = -uy, ux

                for distance in distances:
                    waypoint = {"x": mx + nx * distance, "y": my + ny * distance}
                    override[field] = waypoint
                    candidate = _all_compiled_corridors(model, channel_spacing)
                    quality = _route_defect_score(model, candidate, obstacle_clearance)
                    objective = quality + (abs(distance), reaction.id, edge_key,
                                 field, distance)
                    if best is None or objective < best[0]:
                        best = (objective, reaction, edge_key, field,
                                dict(waypoint), distance)
                override.pop(field, None)

        if best is None or best[0][:3] >= current_quality:
            break
        _, reaction, edge_key, field, waypoint, distance = best
        override = _map_info(reaction).setdefault(
            "edge_route_overrides", {}
        ).setdefault(edge_key, {})
        override[field] = waypoint
        changes.append({
            "reaction": reaction.id,
            "edge": edge_key,
            "field": field,
            "waypoint": waypoint,
            "detour": distance,
        })

    final = _all_compiled_corridors(model, channel_spacing)
    final_quality = _route_defect_score(model, final, obstacle_clearance)
    return {
        "enabled": True,
        "defects_before": before_quality[0],
        "defects_after": final_quality[0],
        "crossings_before": before_quality[1],
        "crossings_after": final_quality[1],
        "obstacles_before": before_quality[2],
        "obstacles_after": final_quality[2],
        "obstacle_hits_after": _route_obstacle_hits(
            model, final, obstacle_clearance
        ),
        "changes": changes,
    }


def _optimize_flexible_nodes(model, channel_spacing, max_nudge=160, step=40,
                             max_iterations=6, obstacle_clearance=42.0):
    """Move explicitly flexible metabolite nodes only when route defects drop."""
    before_edges = _all_compiled_corridors(model, channel_spacing)
    before_quality = _route_defect_score(model, before_edges, obstacle_clearance)
    changes = []

    offsets = [(dx, dy)
               for dx in range(-max_nudge, max_nudge + 1, step)
               for dy in range(-max_nudge, max_nudge + 1, step)
               if dx or dy]
    offsets.sort(key=lambda value: (math.hypot(*value), abs(value[0]) + abs(value[1]),
                                    value[0], value[1]))

    for _ in range(max_iterations):
        compiled = _all_compiled_corridors(model, channel_spacing)
        current_quality = _route_defect_score(model, compiled, obstacle_clearance)
        obstacle_ids = {
            obstacle
            for item in _route_obstacle_hits(model, compiled, obstacle_clearance)
            for obstacle in item["obstacles"]
        }
        candidates = [
            metabolite for metabolite in model.metabolites
            if metabolite.id in obstacle_ids
            and _map_info(metabolite).get("route_node_flexible")
        ]
        if not candidates:
            break

        best = None
        for metabolite in sorted(candidates, key=lambda item: item.id):
            info = _map_info(metabolite)
            origin_x, origin_y = float(info["x"]), float(info["y"])
            for dx, dy in offsets:
                info["x"], info["y"] = origin_x + dx, origin_y + dy
                candidate = _all_compiled_corridors(model, channel_spacing)
                quality = _route_defect_score(model, candidate, obstacle_clearance)
                objective = quality + (math.hypot(dx, dy), abs(dx) + abs(dy),
                                       metabolite.id, dx, dy)
                if best is None or objective < best[0]:
                    best = (objective, metabolite, origin_x, origin_y, dx, dy)
            info["x"], info["y"] = origin_x, origin_y

        if best is None or best[0][:3] >= current_quality:
            break
        _, metabolite, origin_x, origin_y, dx, dy = best
        info = _map_info(metabolite)
        info["x"], info["y"] = origin_x + dx, origin_y + dy
        changes.append({
            "metabolite": metabolite.id,
            "from": [origin_x, origin_y],
            "to": [info["x"], info["y"]],
        })

    final = _all_compiled_corridors(model, channel_spacing)
    final_quality = _route_defect_score(model, final, obstacle_clearance)
    return {
        "enabled": True,
        "defects_before": before_quality[0],
        "defects_after": final_quality[0],
        "crossings_before": before_quality[1],
        "crossings_after": final_quality[1],
        "obstacles_before": before_quality[2],
        "obstacles_after": final_quality[2],
        "changes": changes,
    }


def _compiled_corridors_for_reaction(reaction, channel_spacing, channel_override=None):
    """Compile the lightweight corridor approximation used by route QA.

    This mirrors the per-edge channel semantics written into edge_routes
    without mutating model state. It lets the channel optimizer evaluate
    candidate assignments before committing them.
    """
    info = _map_info(reaction)
    if info.get("hidden") or not {"x", "y"} <= info.keys():
        return []
    style = info.get("route_style", "bezier")
    if style == "reaction_bus":
        bus = _reaction_bus_geometry(
            reaction,
            float(info.get("route_bus_half_length", 34.0)),
        )
        if bus is None:
            return []
        compiled = []
        for source, target in _edge_pairs(reaction):
            source_info, target_info = _map_info(source), _map_info(target)
            compiled.append({
                "reaction": reaction.id,
                "source": source.id,
                "target": target.id,
                "points": [
                    (float(source_info["x"]), float(source_info["y"])),
                    (float(bus["input"][0]), float(bus["input"][1])),
                    (float(bus["output"][0]), float(bus["output"][1])),
                    (float(target_info["x"]), float(target_info["y"])),
                ],
            })
        return compiled

    channel = int(info.get("route_channel", 0) if channel_override is None else channel_override)
    spacing = float(info.get("route_channel_spacing", channel_spacing))
    default_offset = channel * spacing if style == "bundle_spline" else 0.0
    overrides = info.get("edge_route_overrides") or {}
    rx, ry = float(info["x"]), float(info["y"])
    compiled = []
    for source, target in _edge_pairs(reaction):
        source_info, target_info = _map_info(source), _map_info(target)
        sx, sy = float(source_info["x"]), float(source_info["y"])
        tx, ty = float(target_info["x"]), float(target_info["y"])
        key = f"{source.id}->{target.id}"
        override = overrides.get(key, {})
        edge_style = override.get("style", style)
        if edge_style == "arc":
            edge_style = "bundle_spline"
        edge_channel = int(override.get("channel", channel))
        edge_offset = float(override.get(
            "channel_offset",
            edge_channel * spacing if edge_style == "bundle_spline" else default_offset,
        ))
        ux, uy = _unit(tx - sx, ty - sy)
        hx, hy = rx - uy * edge_offset, ry + ux * edge_offset
        points = [(sx, sy)]
        before = override.get("waypoint_before")
        after = override.get("waypoint_after")
        if before:
            points.append((float(before["x"]), float(before["y"])))
        points.append((hx, hy))
        if after:
            points.append((float(after["x"]), float(after["y"])))
        points.append((tx, ty))
        compiled.append({
            "reaction": reaction.id,
            "source": source.id,
            "target": target.id,
            "points": points,
        })
    return compiled


def _symmetric_channels(count):
    order = [0]
    for index in range(1, count + 2):
        order.extend([-index, index])
    return order


def _optimize_bundle_channels(model, bundle_members, channel_spacing, max_exact_size=7):
    """Choose bundle channels by minimizing observable corridor crossings.

    Explicit channels are starting hints unless route_channel_locked is set.
    Small bundles use exact search; larger bundles retain their deterministic
    current assignment instead of introducing an opaque heuristic.
    """
    results = {}
    for bundle, reactions in sorted(bundle_members.items()):
        mutable = [reaction for reaction in reactions
                   if not _map_info(reaction).get("route_channel_locked")]
        if len(mutable) <= 1 or len(mutable) > max_exact_size:
            results[bundle] = {
                "optimized": False,
                "reason": "bundle size outside exact-search range",
            }
            continue

        mutable_ids = {reaction.id for reaction in mutable}
        locked_channels = {
            int(_map_info(reaction).get("route_channel", 0))
            for reaction in reactions if reaction.id not in mutable_ids
        }
        before_assignment = {
            reaction.id: int(_map_info(reaction).get("route_channel", 0))
            for reaction in mutable
        }

        pool = []
        for reaction in mutable:
            value = before_assignment[reaction.id]
            if value not in locked_channels and value not in pool:
                pool.append(value)
        for value in _symmetric_channels(len(reactions) + len(locked_channels)):
            if value not in locked_channels and value not in pool:
                pool.append(value)
            if len(pool) >= len(mutable):
                break
        pool = pool[:len(mutable)]

        fixed_edges = []
        for reaction in model.reactions:
            if reaction.id not in mutable_ids:
                fixed_edges.extend(_compiled_corridors_for_reaction(
                    reaction, channel_spacing
                ))

        current_edges = list(fixed_edges)
        for reaction in mutable:
            current_edges.extend(_compiled_corridors_for_reaction(
                reaction, channel_spacing, before_assignment[reaction.id]
            ))
        before_crossings = _estimate_route_crossings(current_edges)

        best = None
        for permutation in itertools.permutations(pool, len(mutable)):
            assignment = dict(zip((reaction.id for reaction in mutable), permutation))
            candidate_edges = list(fixed_edges)
            for reaction in mutable:
                candidate_edges.extend(_compiled_corridors_for_reaction(
                    reaction, channel_spacing, assignment[reaction.id]
                ))
            crossings = _estimate_route_crossings(candidate_edges)
            movement = sum(
                abs(assignment[reaction.id] - before_assignment[reaction.id])
                for reaction in mutable
            )
            centeredness = sum(abs(value) for value in permutation)
            score = (crossings, movement, centeredness, permutation)
            if best is None or score < best[0]:
                best = (score, assignment, candidate_edges)

        _, assignment, best_edges = best
        for reaction in mutable:
            _map_info(reaction)["route_channel"] = int(assignment[reaction.id])
        results[bundle] = {
            "optimized": True,
            "crossings_before": before_crossings,
            "crossings_after": _estimate_route_crossings(best_edges),
            "channels_before": before_assignment,
            "channels_after": assignment,
        }
    return results


def _compile_route_grammar(model, policy):
    """Compile semantic routing intent into per-edge route specifications.

    The compiler adds three missing concepts to the generic D3 renderer:

    * node ports (N/S/E/W) so curves leave/enter nodes coherently;
    * bundle channels so nearby long edges do not occupy the same corridor;
    * cubic handle policy so non-main edges are smooth rather than rigid.

    Manual `route_channel` remains available as an explicit scientific-layout
    constraint.  Otherwise channels are assigned deterministically within each
    bundle around zero.
    """
    enabled = policy.get("enabled", True)
    if not enabled:
        return {"enabled": False, "edge_routes": 0, "bundles": {}, "estimated_crossings": 0}

    channel_spacing = float(policy.get("channel_spacing", 22.0))
    handle_scale = float(policy.get("handle_scale", 0.30))
    obstacle_clearance = float(policy.get("obstacle_clearance", 42.0))
    multi_metabolite_style = str(policy.get("multi_metabolite_style", "paired"))
    reaction_bus_half_length = float(policy.get("reaction_bus_half_length", 34.0))
    if channel_spacing <= 0 or not (0.05 <= handle_scale <= 0.8):
        raise ValueError("route_policy channel_spacing/handle_scale out of range")
    if multi_metabolite_style not in {"paired", "reaction_bus"}:
        raise ValueError("route_policy multi_metabolite_style must be paired/reaction_bus")
    if reaction_bus_half_length <= 0:
        raise ValueError("route_policy reaction_bus_half_length must be positive")

    flexible_origins = {
        metabolite.id: (
            float(_map_info(metabolite)["x"]),
            float(_map_info(metabolite)["y"]),
        )
        for metabolite in model.metabolites
        if _map_info(metabolite).get("route_node_flexible")
        and {"x", "y"} <= _map_info(metabolite).keys()
    }

    bundle_members = {}
    for reaction in model.reactions:
        info = _map_info(reaction)
        if info.get("hidden"):
            continue
        role = info.get("scientific_role", "other")
        style = info.get("route_style")
        reactants, products = _visible_reaction_sides(reaction)
        is_multi_endpoint = len(reactants) > 1 or len(products) > 1
        if (multi_metabolite_style == "reaction_bus"
                and is_multi_endpoint
                and not info.get("route_style_locked")):
            style = "reaction_bus"
        elif style == "reaction_bus":
            style = "reaction_bus"
        elif style in {"arc", "bundle_spline"}:
            style = "bundle_spline"
        elif style in {"cubic", "cubic_spline"}:
            style = "cubic_spline"
        elif not style:
            style = "straight" if role in {"main"} else "cubic_spline"
        info["route_style"] = style
        if style == "reaction_bus":
            info["route_bus_half_length"] = reaction_bus_half_length
        if style == "bundle_spline":
            bundle = info.get("route_bundle") or role
            bundle_members.setdefault(bundle, []).append(reaction)

    bundle_summary = {}
    for bundle, reactions in sorted(bundle_members.items()):
        explicit = [r for r in reactions if "route_channel" in _map_info(r)]
        implicit = [r for r in reactions if "route_channel" not in _map_info(r)]
        occupied = {int(_map_info(r)["route_channel"]) for r in explicit}
        # Symmetric channel order keeps the whole bundle centered.
        order = _symmetric_channels(len(reactions))
        available = [value for value in order if value not in occupied]
        for reaction, channel in zip(sorted(implicit, key=lambda r: (_map_info(r).get("y", 0), r.id)), available):
            _map_info(reaction)["route_channel"] = channel

    optimization = {}
    if policy.get("optimize_channels", False):
        optimization = _optimize_bundle_channels(
            model,
            bundle_members,
            channel_spacing,
            max_exact_size=int(policy.get("max_exact_bundle_size", 7)),
        )

    for bundle, reactions in sorted(bundle_members.items()):
        bundle_summary[bundle] = {
            "reactions": len(reactions),
            "channels": {r.id: int(_map_info(r).get("route_channel", 0)) for r in reactions},
        }
        if bundle in optimization:
            bundle_summary[bundle]["optimization"] = optimization[bundle]

    edge_optimization = {"enabled": False, "crossings_before": None,
                         "crossings_after": None, "changes": []}
    if policy.get("optimize_edge_routes", False):
        edge_optimization = _optimize_edge_routes(
            model,
            channel_spacing,
            max_channel=int(policy.get("max_edge_channel", 6)),
            max_iterations=int(policy.get("max_edge_iterations", 12)),
            obstacle_clearance=obstacle_clearance,
        )

    hub_optimization = {"enabled": False, "crossings_before": None,
                        "crossings_after": None, "changes": []}
    if policy.get("optimize_route_hubs", False):
        hub_optimization = _optimize_route_hubs(
            model,
            channel_spacing,
            max_nudge=int(policy.get("max_hub_nudge", 160)),
            step=int(policy.get("hub_nudge_step", 40)),
            max_iterations=int(policy.get("max_hub_iterations", 6)),
            obstacle_clearance=obstacle_clearance,
        )

    waypoint_optimization = {"enabled": False, "crossings_before": None,
                             "crossings_after": None, "changes": []}
    if policy.get("optimize_route_waypoints", False):
        waypoint_optimization = _optimize_route_waypoints(
            model,
            channel_spacing,
            max_detour=int(policy.get("max_waypoint_detour", 240)),
            detour_step=int(policy.get("waypoint_detour_step", 40)),
            max_iterations=int(policy.get("max_waypoint_iterations", 8)),
            obstacle_clearance=obstacle_clearance,
        )

    node_optimization = {"enabled": False, "defects_before": None,
                         "defects_after": None, "changes": []}
    if policy.get("optimize_flexible_nodes", False):
        node_optimization = _optimize_flexible_nodes(
            model,
            channel_spacing,
            max_nudge=int(policy.get("max_node_nudge", 160)),
            step=int(policy.get("node_nudge_step", 40)),
            max_iterations=int(policy.get("max_node_iterations", 6)),
            obstacle_clearance=obstacle_clearance,
        )

    compaction = {
        "enabled": False,
        "applied": False,
        "changes": [],
    }
    if policy.get("compact_zero_defect_routes", False):
        compaction = _compact_zero_defect_routes(
            model,
            channel_spacing,
            flexible_origins,
            obstacle_clearance=obstacle_clearance,
            max_node_nudge=int(policy.get("max_node_nudge", 160)),
            node_step=int(policy.get("node_nudge_step", 40)),
            max_edge_channel=int(policy.get("max_edge_channel", 6)),
            max_passes=int(policy.get("max_compaction_passes", 4)),
        )

    compiled = []
    route_count = 0
    reaction_bus_count = 0
    for reaction in model.reactions:
        info = _map_info(reaction)
        if info.get("hidden") or not {"x", "y"} <= info.keys():
            continue
        style = info.get("route_style", "bezier")
        bus = None
        if style == "reaction_bus":
            bus = _reaction_bus_geometry(
                reaction,
                float(info.get("route_bus_half_length", reaction_bus_half_length)),
            )
            if bus is not None:
                info["reaction_bus"] = bus
                reaction_bus_count += 1
        channel = int(info.get("route_channel", 0))
        spacing = float(info.get("route_channel_spacing", channel_spacing))
        default_offset = channel * spacing if style == "bundle_spline" else 0.0
        edge_routes = {}
        overrides = info.get("edge_route_overrides") or {}
        for source, target in _edge_pairs(reaction):
            source_info, target_info = _map_info(source), _map_info(target)
            sx, sy = float(source_info["x"]), float(source_info["y"])
            tx, ty = float(target_info["x"]), float(target_info["y"])
            rx, ry = float(info["x"]), float(info["y"])
            key = f"{source.id}->{target.id}"
            override = overrides.get(key, {})
            if bus is not None:
                input_x, input_y = map(float, bus["input"])
                output_x, output_y = map(float, bus["output"])
                source_port = _cardinal_port(input_x - sx, input_y - sy)
                target_port = _cardinal_port(output_x - tx, output_y - ty)
                source_vector = _unit(input_x - sx, input_y - sy)
                target_vector = _unit(output_x - tx, output_y - ty)
                # Hyperedge coherence is an invariant: per-pair overrides may
                # not split a shared reaction bus into independent corridors.
                edge_style = "reaction_bus"
                edge_channel = 0
                edge_offset = 0.0
            else:
                source_port = _cardinal_port(rx - sx, ry - sy)
                target_port = _cardinal_port(rx - tx, ry - ty)
                source_vector = _unit(rx - sx, ry - sy)
                target_vector = _unit(rx - tx, ry - ty)
                edge_style = override.get("style", style)
                if edge_style == "arc":
                    edge_style = "bundle_spline"
                edge_channel = int(override.get("channel", channel))
                edge_offset = float(override.get(
                    "channel_offset",
                    edge_channel * spacing if edge_style == "bundle_spline" else default_offset,
                ))
            edge = {
                "style": edge_style,
                "source_port": override.get("source_port", source_port),
                "target_port": override.get("target_port", target_port),
                "source_vector": list(override.get("source_vector", source_vector)),
                "target_vector": list(override.get("target_vector", target_vector)),
                "channel": edge_channel,
                "channel_offset": edge_offset,
                "handle_scale": float(override.get("handle_scale", info.get("route_handle_scale", handle_scale))),
            }
            if bus is None and override.get("waypoint_before"):
                edge["waypoint_before"] = dict(override["waypoint_before"])
            if bus is None and override.get("waypoint_after"):
                edge["waypoint_after"] = dict(override["waypoint_after"])
            edge_routes[key] = edge
            route_count += 1

            points = [(sx, sy)]
            if bus is not None:
                points.append((float(bus["input"][0]), float(bus["input"][1])))
                points.append((float(bus["output"][0]), float(bus["output"][1])))
            else:
                # Approximate the rendered corridor by shifting the hidden
                # reaction hub perpendicular to the endpoint chord.
                ux, uy = _unit(tx - sx, ty - sy)
                hx, hy = rx - uy * edge_offset, ry + ux * edge_offset
                if edge.get("waypoint_before"):
                    points.append((float(edge["waypoint_before"]["x"]),
                                   float(edge["waypoint_before"]["y"])))
                points.append((hx, hy))
                if edge.get("waypoint_after"):
                    points.append((float(edge["waypoint_after"]["x"]),
                                   float(edge["waypoint_after"]["y"])))
            points.append((tx, ty))
            compiled.append({
                "reaction": reaction.id,
                "source": source.id,
                "target": target.id,
                "points": points,
            })
        if edge_routes:
            info["edge_routes"] = edge_routes

    crossing_pairs = _route_crossing_pairs(compiled)
    obstacle_hits = _route_obstacle_hits(
        model, compiled, obstacle_clearance
    )
    return {
        "enabled": True,
        "edge_routes": route_count,
        "reaction_buses": reaction_bus_count,
        "bundles": bundle_summary,
        "edge_optimization": edge_optimization,
        "hub_optimization": hub_optimization,
        "waypoint_optimization": waypoint_optimization,
        "node_optimization": node_optimization,
        "compaction": compaction,
        "channel_spacing": channel_spacing,
        "handle_scale": handle_scale,
        "estimated_crossings": len(crossing_pairs),
        "crossing_pairs": crossing_pairs,
        "estimated_obstacle_hits": len(obstacle_hits),
        "obstacle_hits": obstacle_hits,
    }


def _estimate_label_box(obj, info, x, y, kind, font_size):
    name = str(info.get("display_name") or obj.id)
    lines = name.split("\n")
    if kind == "reactions":
        lines.append("0.000")  # reserve one flux-value line
    width = max((len(line) for line in lines), default=1) * font_size * 0.56
    height = len(lines) * font_size * 1.25
    dx = float(info.get("label_dx", 0.0))
    dy = float(info.get("label_dy", 0.0))
    cx, cy = x + dx, y + dy
    return (cx - width / 2, cy - height / 2, cx + width / 2, cy + height / 2)


def _overlap(a, b, padding=4.0):
    return not (a[2] + padding < b[0] or b[2] + padding < a[0]
                or a[3] + padding < b[1] or b[3] + padding < a[1])


def _resolve_label_offsets(model, policy):
    """Small deterministic collision reducer for publication labels.

    It is intentionally conservative: explicit label offsets win; otherwise it
    tries a short, deterministic set of offsets and records the remaining
    approximate collision count for QA.
    """
    if not policy.get("enabled", True):
        return {"auto_placed": 0, "approx_collisions": 0}
    met_font = float(policy.get("metabolite_font", 15))
    rxn_font = float(policy.get("reaction_font", 9))
    candidates = policy.get("candidate_offsets") or [
        [0, -18], [0, 18], [18, -12], [-18, -12], [26, 8], [-26, 8], [0, -34], [0, 34]
    ]
    placed = []
    auto = 0
    ordered = []
    for kind in ("metabolites", "reactions"):
        for obj in getattr(model, kind):
            info = _map_info(obj)
            if info.get("hidden") or not {"x", "y"} <= info.keys():
                continue
            priority = int(info.get("label_priority", 50 if kind == "metabolites" else 20))
            ordered.append((-priority, kind, obj))
    ordered.sort(key=lambda item: (item[0], item[1], item[2].id))

    for _, kind, obj in ordered:
        info = _map_info(obj)
        font = met_font if kind == "metabolites" else rxn_font
        if "label_dx" in info or "label_dy" in info:
            placed.append((obj.id, _estimate_label_box(obj, info, info["x"], info["y"], kind, font)))
            continue
        chosen = candidates[0]
        for dx, dy in candidates:
            info["label_dx"], info["label_dy"] = float(dx), float(dy)
            box = _estimate_label_box(obj, info, info["x"], info["y"], kind, font)
            if not any(_overlap(box, other) for _, other in placed):
                chosen = [dx, dy]
                break
        info["label_dx"], info["label_dy"] = map(float, chosen)
        placed.append((obj.id, _estimate_label_box(obj, info, info["x"], info["y"], kind, font)))
        auto += 1

    collisions = 0
    for i, (_, a) in enumerate(placed):
        for _, b in placed[i + 1:]:
            collisions += int(_overlap(a, b))
    return {"auto_placed": auto, "approx_collisions": collisions}


def _build_structure_checks(model, structures, position_tolerance=2.0):
    lookup = _lookup(model)
    checks = []
    for index, structure in enumerate(structures or []):
        name = structure.get("name") or f"structure_{index + 1}"
        center = structure.get("center") or [0.0, 0.0]
        cx, cy = map(float, center)
        rx, ry = float(structure.get("rx", 1.0)), float(structure.get("ry", 1.0))
        errors = []
        angles = []
        missing = []
        ordered_items = _structure_items(structure)
        ordered_items.sort(key=lambda item: _map_info(lookup[item["id"]][1]).get(
            "layout_structure_order_index", 10**9
        ) if item["id"] in lookup else 10**9)
        for item in ordered_items:
            item_id = item["id"]
            if item_id not in lookup:
                missing.append(item_id)
                continue
            info = _map_info(lookup[item_id][1])
            if not {"x", "y", "layout_theta_deg"} <= info.keys():
                missing.append(item_id)
                continue
            angle = float(info["layout_theta_deg"])
            theta = math.radians(angle)
            base_dx, base_dy = rx * math.cos(theta), ry * math.sin(theta)
            radial_length = math.hypot(base_dx, base_dy) or 1.0
            radial_offset = float(info.get("layout_radial_offset", 0.0))
            expected = (
                cx + base_dx + base_dx / radial_length * radial_offset,
                cy + base_dy + base_dy / radial_length * radial_offset,
            )
            errors.append(math.hypot(float(info["x"]) - expected[0],
                                     float(info["y"]) - expected[1]))
            angles.append(angle)
        if len(angles) < 2:
            monotonic = True
        else:
            increasing = angles[-1] >= angles[0]
            monotonic = all(
                (b >= a - 1e-9) if increasing else (b <= a + 1e-9)
                for a, b in zip(angles, angles[1:])
            )
        max_error = max(errors, default=0.0)
        checks.append({
            "name": name,
            "type": structure.get("type", "ellipse"),
            "items": len(angles),
            "missing": missing,
            "angular_monotonic": monotonic,
            "max_position_error": max_error,
            "position_tolerance": position_tolerance,
            "pass": not missing and monotonic and max_error <= position_tolerance,
        })
    return checks


def build_publication_qa(model, layout):
    """Return cheap, deterministic checks for the final visual scene."""
    if layout.get('circular_backbone', {}).get('enabled'):
        from circular_layout import circular_publication_qa
        return circular_publication_qa(model, layout)
    width, height = map(float, layout["canvas"])
    lookup = _lookup(model)
    missing_coordinates = []
    out_of_bounds = []
    visible = 0
    role_counts = Counter()
    for item_id, (kind, obj) in lookup.items():
        info = _map_info(obj)
        if info.get("hidden"):
            continue
        visible += 1
        if not {"x", "y"} <= info.keys() or not all(math.isfinite(float(info[k])) for k in ("x", "y")):
            missing_coordinates.append(item_id)
            continue
        x, y = float(info["x"]), float(info["y"])
        if not (0 <= x <= width and 0 <= y <= height):
            out_of_bounds.append(item_id)
        if info.get("scientific_role"):
            role_counts[info["scientific_role"]] += 1

    lane_checks = []
    spec = layout.get("scientific_layout") or {}
    for lane in spec.get("lanes", []):
        ids = lane.get("items") or []
        axis = lane.get("axis", "x")
        coords = []
        for item_id in ids:
            if item_id in lookup:
                info = _map_info(lookup[item_id][1])
                coords.append(float(info.get(axis, float("nan"))))
        monotonic = all(math.isfinite(v) for v in coords) and all(a < b for a, b in zip(coords, coords[1:]))
        lane_checks.append({"name": lane.get("name", "unnamed"), "axis": axis,
                            "strictly_monotonic": monotonic, "items": len(coords)})

    qa_policy = layout.get("publication_qa_policy") or {}
    structure_checks = _build_structure_checks(
        model,
        spec.get("structures", []),
        position_tolerance=float(qa_policy.get("structure_position_tolerance", 2.0)),
    )
    route_policy = spec.get("route_policy") or {}
    route_visual = _route_visual_metrics(
        model,
        float(route_policy.get("channel_spacing", 22.0)),
        edge_clearance=float(qa_policy.get("edge_edge_clearance", 18.0)),
        max_route_stretch=float(qa_policy.get("max_route_stretch", 2.4)),
    )
    enforce_visual = bool(qa_policy.get("enforce_visual_metrics", False))
    visual_pass = (route_visual["near_edge_count"] == 0
                   and route_visual["route_stretch_violation_count"] == 0)

    return {
        "visible_objects": visible,
        "missing_coordinates": missing_coordinates,
        "out_of_bounds": out_of_bounds,
        "role_counts": dict(sorted(role_counts.items())),
        "lane_checks": lane_checks,
        "structure_checks": structure_checks,
        "route_visual": route_visual,
        "pass": not missing_coordinates and not out_of_bounds
                and all(check["strictly_monotonic"] for check in lane_checks)
                and all(check["pass"] for check in structure_checks)
                and (visual_pass or not enforce_visual),
    }
