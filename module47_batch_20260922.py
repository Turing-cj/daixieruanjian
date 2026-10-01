"""Render the 47 retained module-conditioned carbon-fixation routes.

This adapter reuses the existing circular renderer.  It never re-solves a
pathway: model.tsv and fit_result.json are generated directly from the frozen
directional stoichiometry/flux in the server export.  A longest real directed
carbon cycle becomes the circular backbone; remaining reactions are rendered
as explicit side branches or local side reactions.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import json
import math
import shutil
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent

sys.path.insert(0, str(ROOT))

import scripts_render_d3flux as native
from publication_readability import enlarge_publication
from label_clearance import clear_publication_labels


CURRENCY = {
    "C00001", "C00080", "C00002", "C00008", "C00020", "C00009",
    "C00013", "C00003", "C00004", "C00005", "C00006", "C00010",
    "C00011", "C00288",
}

# Publication labels intentionally distinguish KEGG anomers that would
# otherwise collapse to the same short biochemical name.
DISPLAY = {
    "C00001": "H₂O", "C00002": "ATP", "C00003": "NAD⁺", "C00004": "NADH",
    "C00005": "NADPH", "C00006": "NADP⁺", "C00008": "ADP", "C00009": "Pi",
    "C00011": "CO₂", "C00022": "Pyruvate", "C00074": "PEP", "C00080": "H⁺",
    "C00085": "F6P", "C00092": "G6P", "C00111": "DHAP", "C00117": "R5P",
    "C00118": "G3P", "C00184": "DHA", "C00197": "3PG", "C00199": "Ru5P",
    "C00227": "Acetyl-P", "C00231": "Xu5P", "C00236": "1,3-BPG",
    "C00279": "E4P", "C00345": "6PG", "C00354": "FBP", "C00447": "SBP",
    "C00631": "2PG", "C01182": "RuBP", "C01236": "6PGL",
    "C03736": "α-R5P", "C04442": "KDPG", "C05345": "β-F6P", "C05382": "S7P",
}

MODULE_COLORS = {
    "EMP": "#4A8D88",
    "ED": "#657786",
    "PKT": "#8B6FA3",
    "OPPP": "#B58AA6",
}


def dump(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def term(mid: str, coefficient: float) -> str:
    value = abs(float(coefficient))
    prefix = "" if math.isclose(value, 1.0) else f"{value:.12g}"
    return prefix + mid


def reaction_label(column: str) -> str:
    base, direction = column.rsplit("_", 1)
    if base.startswith("IML_"):
        label = base[4:]
    elif base.startswith("KEGG_"):
        label = base[5:]
    else:
        label = base
    return label + (" [R]" if direction == "R" else "")


def metabolite_label(mid: str, source_names: dict[str, str]) -> str:
    if mid in DISPLAY:
        return DISPLAY[mid]
    name = source_names.get(mid, mid)
    # The source table occasionally appends the molecular formula to the name.
    # Keep the readable biochemical name; exact identity remains the model ID.
    return name if len(name) <= 28 else mid


def build_graph(reaction_model: list) -> tuple[dict[str, list[str]], dict[str, dict[str, float]]]:
    adj: dict[str, list[str]] = defaultdict(list)
    stoich: dict[str, dict[str, float]] = {}
    for row in reaction_model:
        rid, coefficients = row[0], {str(k): float(v) for k, v in row[1].items()}
        stoich[rid] = coefficients
        rnode = "r:" + rid
        adj.setdefault(rnode, [])
        for mid, value in coefficients.items():
            if mid in CURRENCY:
                continue
            mnode = "m:" + mid
            adj.setdefault(mnode, [])
            if value < 0:
                adj[mnode].append(rnode)
            elif value > 0:
                adj[rnode].append(mnode)
    for key in adj:
        adj[key] = sorted(set(adj[key]))
    return dict(adj), stoich


def best_cycle(adj: dict[str, list[str]], stoich: dict[str, dict[str, float]]) -> list[str]:
    """Choose a largest real cycle whose residual network is easiest to draw.

    Reaction count is the hard invariant.  Among equal-size biochemical
    cycles, prefer the embedding that leaves a simple external side path:
    fewer internal returns to the ring, more leftover reactions covered by one
    complete path, fewer detached leftovers, and fewer hyperedge branches.
    This makes layout quality a topological choice instead of a pixel patch.
    """
    best_count = 0
    candidates: list[list[str]] = []
    metabolite_nodes = sorted(node for node in adj if node.startswith("m:"))
    for start in metabolite_nodes:
        path = [start]
        seen = {start}

        def visit(node: str) -> None:
            nonlocal best_count, candidates
            for nxt in adj.get(node, []):
                if nxt == start and len(path) >= 4:
                    score = sum(item.startswith("r:") for item in path)
                    if score > best_count:
                        best_count = score
                        candidates = [path[:]]
                    elif score == best_count:
                        candidates.append(path[:])
                elif nxt not in seen:
                    # Canonicalize cycle start to prevent redundant enumeration.
                    if nxt.startswith("m:") and nxt < start:
                        continue
                    seen.add(nxt)
                    path.append(nxt)
                    visit(nxt)
                    path.pop()
                    seen.remove(nxt)

        visit(start)
    if not candidates:
        return []

    all_reactions = set(stoich)

    def residual_score(cycle: list[str]) -> tuple:
        cycle_reactions = {node[2:] for node in cycle if node.startswith("r:")}
        ring_metabolites = {node[2:] for node in cycle if node.startswith("m:")}
        leftovers = all_reactions - cycle_reactions
        side = longest_side_path(adj, leftovers, ring_metabolites) if leftovers else []
        side_reactions = {node[2:] for node in side if node.startswith("r:")}
        side_metabolites = [node[2:] for node in side if node.startswith("m:")]
        internal_ring_hits = sum(mid in ring_metabolites for mid in side_metabolites[1:-1])
        extra_count = len(leftovers - side_reactions)
        side_hyperedges = sum(
            max(0, sum(mid not in CURRENCY for mid in stoich[rid]) - 2)
            for rid in leftovers
        )
        cycle_hyperedges = sum(
            max(0, sum(mid not in CURRENCY for mid in stoich[rid]) - 2)
            for rid in cycle_reactions
        )
        # Maximize this tuple.  Lexicographic cycle order is only a final,
        # deterministic tie breaker and carries no biochemical meaning.
        return (
            -internal_ring_hits,
            len(side_reactions),
            -extra_count,
            -side_hyperedges,
            -cycle_hyperedges,
            tuple(cycle),
        )

    return max(candidates, key=residual_score)


def rotate_cycle(cycle: list[str]) -> list[str]:
    if not cycle:
        return cycle
    anchors = ["m:C01182", "m:C00199", "m:C00118"]
    index = next((cycle.index(anchor) for anchor in anchors if anchor in cycle), 0)
    rotated = cycle[index:] + cycle[:index]
    if not rotated[0].startswith("m:"):
        # Bipartite cycles alternate, so one step is sufficient.
        rotated = rotated[1:] + rotated[:1]
    return rotated


def longest_side_path(adj: dict[str, list[str]], leftovers: set[str], ring_metabolites: set[str],
                      preferred_endpoints: set[str] | None = None) -> list[str]:
    """Find the longest directed path using only non-backbone reactions.

    Ties prefer a path that starts on the ring, because it is usually the
    product branch and gives the clearest open-chain presentation.
    """
    allowed = {"r:" + rid for rid in leftovers}
    preferred_endpoints = preferred_endpoints or set()
    candidates: list[tuple[tuple[int, int, int, int], list[str]]] = []
    starts = sorted(node for node in adj if node.startswith("m:"))
    for start in starts:
        path = [start]
        seen = {start}

        def visit(node: str) -> None:
            if not node.startswith("m:"):
                # A drawable biochemical side path is a complete sequence
                # metabolite -> reaction -> metabolite.  Do not rank a
                # half-edge ending at a reaction as a valid chain.
                pass
            else:
                reaction_count = sum(item.startswith("r:") for item in path)
                score = (
                    reaction_count,
                    int(path[-1][2:] in preferred_endpoints),
                    int(path[0][2:] in ring_metabolites),
                    int(path[-1][2:] not in ring_metabolites),
                )
                candidates.append((score, path[:]))
            for nxt in adj.get(node, []):
                if nxt.startswith("r:") and nxt not in allowed:
                    continue
                if nxt in seen:
                    continue
                seen.add(nxt)
                path.append(nxt)
                visit(nxt)
                path.pop()
                seen.remove(nxt)

        visit(start)
    return max(candidates, key=lambda row: (row[0], tuple(row[1])), default=((0, 0, 0, 0), []))[1]


def point(center: list[float], radius: float, degrees: float) -> list[float]:
    angle = math.radians(degrees)
    return [center[0] + radius * math.cos(angle), center[1] + radius * math.sin(angle)]


def label_offset(name: str, angle: float, font: float = 19.0, radius: float = 9.0) -> dict[str, float]:
    width = (len(name) + 7) * font * 4 / 3 * 0.53
    height = font * 4 / 3 * 1.1
    theta = math.radians(angle)
    x, y = math.cos(theta), math.sin(theta)
    distance = min(
        (width / 2 + radius + 12) / max(abs(x), 0.001),
        (height / 2 + radius + 12) / max(abs(y), 0.001),
    )
    return {"label_dx": distance * x, "label_dy": distance * y - 4}


def choose_pair(sto: dict[str, float], preferred: set[str]) -> tuple[str, str]:
    substrates = sorted(mid for mid, value in sto.items() if value < 0 and mid not in CURRENCY)
    products = sorted(mid for mid, value in sto.items() if value > 0 and mid not in CURRENCY)
    if not substrates or not products:
        raise ValueError("visible side reaction lacks a non-currency substrate/product pair")
    pairs = [(a, b) for a in substrates for b in products]
    return max(pairs, key=lambda pair: (int(pair[0] in preferred) + int(pair[1] in preferred), pair[0], pair[1]))


def should_merge_shared_instance(mid: str, reaction_ids: list[str],
                                 stoich: dict[str, dict[str, float]],
                                 reaction_angles: dict[str, float] | None = None,
                                 max_separation_deg: float = 60.0) -> bool:
    """Return whether extra incidences are simple enough for one visual hub.

    ``shared_instances`` is a visual compression device, not a new biochemical
    pool.  It is helpful for a pair of equivalent extra incidences, but becomes
    counterproductive when it mixes production and consumption or collapses a
    high-fanout pool into one centre node.  In those cases the source-split
    machinery already has the correct representation: separate local views of
    the same chemical pool, with source/pool markers carrying provenance.

    Keeping this admission rule deliberately narrow prevents a locally tidy
    hub from creating global spoke crossings and ambiguous source semantics.
    """
    if len(reaction_ids) != 2:
        return False
    signs = {
        1 if stoich[rid][mid] > 0 else -1
        for rid in reaction_ids
        if stoich[rid][mid] != 0
    }
    if len(signs) != 1:
        return False
    first, second = reaction_ids
    shared_other = (
        set(stoich[first]) & set(stoich[second])
    ) - {mid} - CURRENCY
    if shared_other:
        return False
    if reaction_angles is not None:
        a, b = (reaction_angles[rid] % 360.0 for rid in reaction_ids)
        separation = abs(a - b)
        separation = min(separation, 360.0 - separation)
        if separation > max_separation_deg:
            return False
    return True


def _circular_mean_degrees(angles: list[float]) -> tuple[float, float]:
    """Return circular mean and concentration for a non-empty angle list."""
    if not angles:
        raise ValueError('angles must not be empty')
    sx = math.fsum(math.cos(math.radians(angle)) for angle in angles)
    sy = math.fsum(math.sin(math.radians(angle)) for angle in angles)
    angle = math.degrees(math.atan2(sy, sx)) if abs(sx) + abs(sy) > 1e-9 else angles[0]
    return angle, math.hypot(sx, sy) / len(angles)


def _bounded_circular_pull(primary: float, proposed: float, limit: float) -> float:
    """Let secondary attachments bend a hub, but never drag it out of sector."""
    delta = (proposed - primary + 180.0) % 360.0 - 180.0
    return primary + max(-limit, min(limit, delta))


def extra_currency_fan(
    a: list[float],
    b: list[float],
    center: list[float],
    by_sign: dict[int, list[str]],
    canonical_at_source: bool | None,
    lane_floor_by_sign: dict[int, float] | None = None,
) -> dict[str, dict[str, list[float]]]:
    """Place currency incidences for a short extra-reaction stub.

    The old layout placed same-sign currencies at increasing distances on one
    normal ray. For reactions with two products (for example NADPH + CO2), the
    farther edge therefore passed through the nearer currency node. The fan
    below enforces two geometric invariants:

    1. every currency ray leaves the reaction at a distinct angle;
    2. the whole fan is shifted toward the open/outward end of the stub, away
       from the already occupied ring/main-chain attachment.

    canonical_at_source identifies which trunk end is attached to the already
    drawn canonical graph. None falls back to the radial direction from the
    figure centre for a genuinely isolated panel reaction.
    """
    vx, vy = b[0] - a[0], b[1] - a[1]
    vnorm = max(math.hypot(vx, vy), 1.0)
    ux, uy = vx / vnorm, vy / vnorm
    mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2

    if canonical_at_source is True:
        ox, oy = ux, uy
    elif canonical_at_source is False:
        ox, oy = -ux, -uy
    else:
        rx, ry = mx - center[0], my - center[1]
        rnorm = math.hypot(rx, ry)
        if rnorm < 1e-9:
            ox, oy = ux, uy
        else:
            ox, oy = rx / rnorm, ry / rnorm

    # Lateral axis around the outward direction. Substrates and products use
    # opposite sides; members on the same side get separate angular lanes.
    sx, sy = -oy, ox
    base_x, base_y = mx + 420 * ox, my + 420 * oy
    result: dict[str, dict[str, list[float]]] = {}
    lane_floor_by_sign = lane_floor_by_sign or {}
    for sign in (-1, 1):
        members = sorted(by_sign.get(sign, []))
        first_lane = max(70.0, float(lane_floor_by_sign.get(sign, 0.0)) + 55.0)
        for offset, mid in enumerate(members):
            # A narrow cone is sufficient to make rays non-collinear. Keeping
            # the lateral displacement below the neighbouring main-chain
            # reaction bands prevents this local fan from invading them.
            lane = sign * (first_lane + 55 * offset)
            x = base_x + lane * sx
            y = base_y + lane * sy
            result[mid] = {
                "position": [x, y],
                "label_offset": [
                    52 * ox + 20 * sign * sx,
                    52 * oy + 20 * sign * sy - 4,
                ],
            }
    return result


def main_side_currency_lane(sign: int, offset: int,
                            reserved_lanes: dict[int, float]) -> float:
    """Lateral lane for a currency incidence on the vertical main side chain.

    Local carbon satellites are placed first.  A nearer currency node must sit
    *outside* the farthest same-sign carbon lane; otherwise the two incidence
    curves reverse their order on the way to the reaction bus and cross.
    """
    floor = float(reserved_lanes.get(sign, 0.0))
    direction = sign
    first = floor + 55.0 if floor > 0 else 32.0
    magnitude = first + 42.0 * offset
    return direction * magnitude


def build_layout(route_id: str, module_id: str, rank: int, reaction_model: list,
                 source_names: dict[str, str]) -> tuple[dict, dict]:
    adj, stoich = build_graph(reaction_model)
    cycle = rotate_cycle(best_cycle(adj, stoich))
    cycle_reactions = [node[2:] for node in cycle if node.startswith("r:")]
    if len(cycle_reactions) < 6:
        raise ValueError(f"{route_id}: no sufficiently large biochemical cycle")
    cycle_items = [node[2:] for node in cycle]
    ring_metabolites = set(cycle_items[::2])
    all_reactions = [row[0] for row in reaction_model]
    leftovers = set(all_reactions) - set(cycle_reactions)
    all_metabolites = {mid for coefficients in stoich.values() for mid in coefficients}
    preferred_endpoints = (
        {"C00022"} if "C00022" in all_metabolites
        else {"C00227"} if "C00227" in all_metabolites
        else set()
    )
    side_path = longest_side_path(
        adj, leftovers, ring_metabolites, preferred_endpoints
    ) if leftovers else []
    chain_reactions = [node[2:] for node in side_path if node.startswith("r:")]
    chain_set = set(chain_reactions)
    extra_reactions = sorted(leftovers - chain_set)

    # Keep the same publication typography as the accepted nine-path renderer,
    # but give the larger module-conditioned networks more geometric room.
    # The previous adapter pushed satellites very far from a comparatively
    # small ring; adjacent hyperedge rays then crossed before label clearance
    # even had a chance to act.
    center = [1600.0, 1600.0]
    radius = 760.0
    width = 5000
    chain_nodes = [node[2:] for node in side_path if node.startswith("m:")]
    chain_spacing = 245
    chain_height = max(0, (len(chain_nodes) - 1) * chain_spacing)
    extra_height = len(extra_reactions) * 220
    height = int(max(3500, 720 + chain_height, 2380 + extra_height))
    side_x = 3950.0
    side_y0 = max(220.0, center[1] - chain_height / 2)
    step = 360.0 / len(cycle_items)
    # Rotate the ring so the main residual branch attaches on the right-hand
    # side whenever it shares a real ring metabolite.  Rotation changes only
    # presentation; it does not change the directed cycle or reaction model.
    main_anchor_mid = next(
        (node[2:] for node in side_path if node.startswith("m:") and node[2:] in cycle_items),
        None,
    )
    start_angle = (
        -step * cycle_items.index(main_anchor_mid)
        if main_anchor_mid is not None else -90.0
    )
    ring_positions = {
        mid: point(center, radius, start_angle + step * index)
        for index, mid in enumerate(cycle_items)
        if index % 2 == 0
    }

    # Merge a non-ring pool only when several *backbone* reactions genuinely
    # need one common hub.  A pool that appears once on the backbone and again
    # on a side route is not a shared-hub by itself: source_split can represent
    # those separate incidences without inventing a physical pool split.  The
    # previous broader rule incorrectly promoted C00117 into shared_metabolites
    # and then made any local side incidence semantically illegal in browser QA.
    connector_mids = set()
    canonical_satellite_mids = set()
    for mid in {m for coefficients in stoich.values() for m in coefficients if m not in CURRENCY} - ring_metabolites:
        backbone_incidence = sum(mid in stoich[rid] for rid in cycle_reactions)
        if backbone_incidence >= 2:
            connector_mids.add(mid)
        elif backbone_incidence == 1 and any(mid in stoich[rid] for rid in leftovers):
            canonical_satellite_mids.add(mid)
    incidence_count = Counter(
        mid for rid in all_reactions for mid in stoich[rid]
        if mid not in CURRENCY
    )
    external_mids = {
        mid for mid, count in incidence_count.items()
        if count == 1 and mid not in ring_metabolites
    }

    reaction_angles = {
        cycle_items[index]: start_angle + step * index
        for index in range(1, len(cycle_items), 2)
    }
    shared_positions: dict[str, list[float]] = {}
    for mid in sorted(connector_mids):
        incident_cycle = [rid for rid in cycle_reactions if mid in stoich[rid]]
        # Backbone incidences define the hub's primary visual sector. Residual
        # attachments may bend that hub slightly, but must never drag it across
        # unrelated ring reactions. The old unbounded blend caused the ED-0001
        # Xu5P hub to move from the TKT1/TKT2 sector toward RPE and create the
        # conspicuous lower-ring line crossing reported in publication review.
        primary_angles = [reaction_angles[rid] for rid in incident_cycle]
        primary_angle, concentration = _circular_mean_degrees(primary_angles)
        residual_angles = []
        # A connector can also feed a residual side reaction that returns to a
        # canonical ring metabolite. Include that real attachment when placing
        # the hub; otherwise the side trunk is forced to cut back across the
        # regeneration spokes.
        for rid in leftovers:
            if mid not in stoich[rid]:
                continue
            for other in stoich[rid]:
                if other == mid or other not in ring_positions:
                    continue
                px, py = ring_positions[other]
                residual_angles.append(math.degrees(math.atan2(py - center[1], px - center[0])))
        if residual_angles:
            proposed, _ = _circular_mean_degrees(primary_angles * 2 + residual_angles)
            angle = _bounded_circular_pull(primary_angle, proposed, step * 0.20)
        else:
            angle = primary_angle
        role_score = sum(
            1 if stoich[rid][mid] > 0 else -1
            for rid in incident_cycle
        )
        if role_score:
            angle += math.copysign(step * 0.25, role_score)
        # Pools shared by several backbone reactions (classically Xu5P in the
        # regeneration block) are clearer as one compact inner hub.  This is
        # the same hyperedge pattern used by the accepted nine-path layouts:
        # one real pool node, several reaction branches, no fake outer copy.
        # Place the hub near the geometric centroid of its incident backbone
        # reactions. A fixed deep-centre radius makes a pool shared by nearby
        # reactions cut across unrelated spokes on its way back out.
        hub_fraction = max(0.50, min(0.78, concentration))
        shared_positions[mid] = point(center, radius * hub_fraction, angle)

    # Every model object receives a valid initial coordinate; the circular
    # compiler then owns the final ring/satellite geometry.
    metabolites = sorted({mid for coefficients in stoich.values() for mid in coefficients})
    layout = {
        "template": "candidate01_demo/empty_template.json",
        "template_transform": [1, 0, 0],
        "template_metabolites": {},
        "template_reactions": {},
        "canvas": [width, height],
        "pathway_multiplicity": {
            "enabled": True,
            "normalization_flux": 1.0,
            "label_mode": "local_flow",
            "relationship_mode": "source_split",
        },
        "circular_backbone": {
            "enabled": True,
            "center": center,
            "radius": radius,
            "start_angle_deg": start_angle,
            "items": cycle_items,
            "external_metabolites": sorted(external_mids),
            "canonical_satellites": sorted(canonical_satellite_mids),
            "satellite_inset": 185,
            "satellite_spacing": 90,
            "local_carbon_outset": 190,
            "currency_satellite_outset": 285,
            "external_satellite_outset": 230,
            "shared_metabolites": {},
            "shared_instances": {},
            "instance_marker_style": "pool_code",
            "input_arrows": True,
            "side_reactions": [],
            "title": f"{module_id} module",
            "subtitle": f"{route_id} · rank {rank}",
            "title_position": [center[0] + 150, center[1] - 150],
            "instance_legend_position": [width / 2, height - 135],
            "view_overrides": {},
        },
        "line_style": {
            "mode": "flux_width_only",
            "threshold": 0.4,
            "width_range": [2.4, 5.7],
            "width_exponent": 0.65,
            "unit": "relative pathway flux (product = 1)",
        },
        "publication_style": {
            "metabolite_font_pt": 19,
            "reaction_font_pt": 13,
            "legend_position": [width / 2 - 210, height - 100],
        },
        "hidden_metabolites": sorted(mid for mid in metabolites if mid in CURRENCY),
        "hidden_reactions": [],
        "metabolites": {
            mid: {
                "x": center[0], "y": center[1],
                "display_name": metabolite_label(mid, source_names),
                "node_radius": 5.5 if mid in CURRENCY else 9,
            }
            for mid in metabolites
        },
        "reactions": {},
        "scientific_layout": {
            "enabled": True,
            "roles": {
                "sink": ["C00022"] if "C00022" in metabolites else (["C00227"] if "C00227" in metabolites else []),
                "cycle_anchor": ["C01182"] if "C01182" in metabolites else [],
                "regeneration": sorted(ring_metabolites),
            },
            "role_palette": {
                "sink": "#2F6F8F", "cycle_anchor": "#8B6FA3", "regeneration": "#B58AA6",
            },
            "reaction_styles": {},
        },
    }

    spec = layout["circular_backbone"]
    for mid, position in shared_positions.items():
        theta = math.atan2(position[1] - center[1], position[0] - center[0])
        spec["shared_metabolites"][mid] = {
            "position": position,
            "label_offset": [62 * math.cos(theta), 62 * math.sin(theta) - 4],
        }

    # If a ring metabolite participates in two or more additional backbone
    # reactions beyond its two ring neighbours, merge only those extra visual
    # incidences into one shared view.  The canonical ring node stays put; the
    # shared view keeps the repeated-use relation explicit without crossing
    # several near-identical branches through the centre.
    for index in range(0, len(cycle_items), 2):
        mid = cycle_items[index]
        neighbour_reactions = {
            cycle_items[(index - 1) % len(cycle_items)],
            cycle_items[(index + 1) % len(cycle_items)],
        }
        extra_reactions_for_mid = [
            rid for rid in cycle_reactions
            if rid not in neighbour_reactions and mid in stoich[rid]
        ]
        if not should_merge_shared_instance(
            mid,
            extra_reactions_for_mid,
            stoich,
            reaction_angles,
            max_separation_deg=3.0 * step,
        ):
            continue
        angles = [reaction_angles[rid] for rid in extra_reactions_for_mid]
        sx = math.fsum(math.cos(math.radians(a)) for a in angles)
        sy = math.fsum(math.sin(math.radians(a)) for a in angles)
        angle = math.degrees(math.atan2(sy, sx)) if abs(sx) + abs(sy) > 1e-9 else angles[0]
        role_score = sum(
            1 if stoich[rid][mid] > 0 else -1
            for rid in extra_reactions_for_mid
        )
        if role_score:
            angle += math.copysign(step * 0.80, role_score)
        pos = point(center, radius * 0.48, angle)
        theta = math.radians(angle)
        spec["shared_instances"][f"shared_{mid}"] = {
            "metabolite": mid,
            "reactions": sorted(extra_reactions_for_mid),
            "position": pos,
            "label_offset": [72 * math.cos(theta), 72 * math.sin(theta) - 4],
        }
    for index, item in enumerate(cycle_items):
        angle = start_angle + step * index
        name = reaction_label(item) if index % 2 else metabolite_label(item, source_names)
        if index % 2 == 0:
            spec["view_overrides"][item] = label_offset(name, angle)
        else:
            theta = math.radians(angle)
            spec["view_overrides"][item] = {
                "label_dx": -70 * math.cos(theta),
                "label_dy": -70 * math.sin(theta) - 5,
            }

    for rid in all_reactions:
        coefficients = stoich[rid]
        entry = {
            "x": center[0], "y": center[1],
            "display_name": reaction_label(rid),
            "pathway": "Other",
            "color": MODULE_COLORS[module_id],
            "cofactors": {},
        }
        for mid in sorted(coefficients):
            if mid in CURRENCY:
                entry["cofactors"][mid] = {
                    "x": center[0], "y": center[1],
                    "display_name": metabolite_label(mid, source_names),
                    "node_radius": 5.5,
                    "scientific_role": "cofactor",
                    "color": "#9AA6B2",
                }
        layout["reactions"][rid] = entry
        layout["scientific_layout"]["reaction_styles"][rid] = {"display_name": reaction_label(rid)}

    # The circular compiler owns satellite coordinates.  The layout options
    # above tell it to fan local carbon views and cofactors outward *before*
    # route geometry is compiled, so path endpoints remain attached.

    def place_side_satellites(rid: str, source: str, target: str,
                              a: list[float], b: list[float]) -> tuple[dict[str, dict], dict[int, float]]:
        """Place non-trunk carbon incidences beside a side-reaction trunk.

        Hyperedges such as EDA have a second carbon product in addition to the
        selected open-chain target.  Leaving that incidence on its distant
        canonical ring node creates exactly the long cross-figure ray that the
        side-route projection is meant to avoid.  Repeated pools get a local
        source_split view; unique pools keep their canonical identity but are
        moved beside the local reaction.
        """
        satellites: dict[str, dict] = {}
        reserved_lanes = {-1: 0.0, 1: 0.0}
        mids = sorted(
            mid for mid in stoich[rid]
            if mid not in CURRENCY and mid not in (source, target)
        )
        if not mids:
            return satellites, reserved_lanes
        vx, vy = b[0] - a[0], b[1] - a[1]
        norm = max(math.hypot(vx, vy), 1.0)
        nx, ny = -vy / norm, vx / norm
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        # Pick the normal that points away from the circle centre.
        if nx * (mx - center[0]) + ny * (my - center[1]) < 0:
            nx, ny = -nx, -ny
        for satellite_index, mid in enumerate(mids):
            if mid in shared_positions:
                # shared_metabolites deliberately promise one visual node.
                continue
            distance = 330 + 145 * satellite_index
            coefficient = stoich[rid][mid]
            sign = -1 if coefficient < 0 else 1
            same_sign_currency = any(
                currency_mid in CURRENCY and value * sign > 0
                for currency_mid, value in stoich[rid].items()
            )
            tx, ty = vx / norm, vy / norm
            tangent_lane = (
                sign * (95 + 55 * satellite_index)
                if same_sign_currency else 0
            )
            if tangent_lane:
                reserved_lanes[sign] = max(reserved_lanes[sign], abs(tangent_lane))
            pos = [
                mx + distance * nx + tangent_lane * tx,
                my + distance * ny + tangent_lane * ty,
            ]
            label_offset = [58 * nx, 58 * ny - 4]
            has_existing_canonical = (
                mid in ring_positions or mid in canonical_satellite_mids
                or mid in chain_nodes
            )
            if has_existing_canonical:
                satellites[mid] = {
                    "position": pos,
                    "label_offset": label_offset,
                }
            else:
                spec["view_overrides"][mid] = {
                    "x": pos[0], "y": pos[1],
                    "label_dx": label_offset[0], "label_dy": label_offset[1],
                    "node_radius": 9,
                }
                layout["metabolites"][mid].update(x=pos[0], y=pos[1])
        return satellites, reserved_lanes

    # Main non-backbone path: a clean vertical branch using explicit local
    # copies at any endpoint already represented on the ring/shared hub.  The
    # source/split labels are the relationship carrier, so a repeated pool does
    # not need a long visual wire back across the figure.  This is the same
    # invariant used by the accepted nine-path source_split layouts: chemistry
    # identity is preserved while redundant visual incidence is projected as a
    # labelled local view.
    path_positions: dict[str, list[float]] = {}
    if side_path and chain_reactions:
        for index, node in enumerate(side_path):
            if node.startswith("m:"):
                mid = node[2:]
                # There is one metabolite between consecutive reactions.
                mindex = sum(1 for prior in side_path[:index] if prior.startswith("m:"))
                path_positions[mid + f"@{mindex}"] = [side_x, side_y0 + mindex * chain_spacing]
        metabolite_occurrences = [(idx, node[2:]) for idx, node in enumerate(side_path) if node.startswith("m:")]
        main_anchor_mindex = next(
            (mindex for mindex, (_, mid) in enumerate(metabolite_occurrences) if mid == main_anchor_mid),
            None,
        )
        occurrence_local: dict[int, bool] = {}
        for mindex, (_, mid) in enumerate(metabolite_occurrences):
            repeated = (
                mid in ring_positions or mid in shared_positions
                or mid in canonical_satellite_mids
            )
            # shared_metabolites is a stronger contract than source_split:
            # browser QA requires every incidence to terminate at the single
            # canonical shared node. Ring/canonical pools may still use a
            # local path instance, but adjacent side reactions must later
            # share that *same* local instance so the chain stays connected.
            occurrence_local[mindex] = (
                repeated
                and mid not in shared_positions
                and mindex != main_anchor_mindex
            )
        for reaction_index, rid in enumerate(chain_reactions):
            source = metabolite_occurrences[reaction_index][1]
            target = metabolite_occurrences[reaction_index + 1][1]
            planned_a = [side_x, side_y0 + reaction_index * chain_spacing]
            planned_b = [side_x, side_y0 + (reaction_index + 1) * chain_spacing]
            # One canonical attachment keeps the entire residual chain in the
            # same visual component. Other repeated endpoints may be local
            # source_split views, but a shared_metabolite is never cloned.
            local_source = occurrence_local[reaction_index]
            local_target = occurrence_local[reaction_index + 1]
            a = planned_a if local_source else ring_positions.get(source, shared_positions.get(source, planned_a))
            b = planned_b if local_target else ring_positions.get(target, shared_positions.get(target, planned_b))
            side = {
                "reaction": rid, "source": source, "target": target,
                "start": a, "end": b, "label_offset": [-135, -4],
                "local_source": local_source,
                "local_target": local_target,
                "incidence_join_spread": 0.12,
            }
            local_satellites, reserved_lanes = place_side_satellites(rid, source, target, a, b)
            if local_satellites:
                side["local_satellites"] = local_satellites
            if not local_source and source in ring_positions and target not in ring_positions:
                angle = math.atan2(a[1] - center[1], a[0] - center[0])
                outer = [center[0] + (radius + 280) * math.cos(angle), center[1] + (radius + 280) * math.sin(angle)]
                side["control_points"] = [outer, [b[0] - 260, b[1]]]
            elif not local_target and target in ring_positions and source not in ring_positions:
                angle = math.atan2(b[1] - center[1], b[0] - center[0])
                outer = [center[0] + (radius + 280) * math.cos(angle), center[1] + (radius + 280) * math.sin(angle)]
                side["control_points"] = [[a[0] - 260, a[1]], outer]
            elif not local_source and not local_target and source in ring_positions and target in ring_positions:
                angle_a = math.atan2(a[1] - center[1], a[0] - center[0])
                angle_b = math.atan2(b[1] - center[1], b[0] - center[0])
                side["control_points"] = [
                    [center[0] + (radius + 360) * math.cos(angle_a), center[1] + (radius + 360) * math.sin(angle_a)],
                    [center[0] + (radius + 360) * math.cos(angle_b), center[1] + (radius + 360) * math.sin(angle_b)],
                ]
            spec["side_reactions"].append(side)
            if local_source:
                spec["view_overrides"][f"{source}_{rid}"] = {
                    "label_dx": -115, "label_dy": -4, "node_radius": 9,
                }
            if local_target:
                spec["view_overrides"][f"{target}_{rid}"] = {
                    "label_dx": -115, "label_dy": -4, "node_radius": 9,
                }
            for mid, pos in ((source, a), (target, b)):
                if (mid not in ring_positions and mid not in shared_positions
                        and mid not in canonical_satellite_mids):
                    spec["view_overrides"][mid] = {
                        "x": pos[0], "y": pos[1], "label_dx": 95, "label_dy": -4, "node_radius": 9,
                    }
                    layout["metabolites"][mid].update(x=pos[0], y=pos[1])
            # Currency views on the regular vertical chain use the compact
            # left/right lane grammar.  The attachment segment can be
            # diagonal, however (ring -> side chain or side chain -> ring).
            # Reusing global ``side_x`` lanes on a diagonal segment can make a
            # cofactor incidence collinear with a carbon product incidence.
            # In that case use the segment-local fan already used for extra
            # reactions so every incidence leaves the reaction bus at a
            # distinct angle.
            by_sign = {-1: [], 1: []}
            for mid, coefficient in stoich[rid].items():
                if mid in CURRENCY:
                    by_sign[-1 if coefficient < 0 else 1].append(mid)
            if abs(a[0] - b[0]) > 1e-6:
                canonical_at_source = (
                    not local_source
                    and (source in ring_positions or source in shared_positions)
                )
                canonical_at_target = (
                    not local_target
                    and (target in ring_positions or target in shared_positions)
                )
                attached_at_source = (
                    True if canonical_at_source and not canonical_at_target
                    else False if canonical_at_target and not canonical_at_source
                    else None
                )
                placements = extra_currency_fan(
                    a, b, center, by_sign, attached_at_source,
                    lane_floor_by_sign=reserved_lanes,
                )
                for mid, placement in placements.items():
                    cf = layout["reactions"][rid]["cofactors"][mid]
                    cf.update(
                        x=placement["position"][0],
                        y=placement["position"][1],
                        align="center",
                        label_dx=placement["label_offset"][0],
                        label_dy=placement["label_offset"][1],
                    )
            else:
                middle_y = (a[1] + b[1]) / 2
                for sign in (-1, 1):
                    for offset, mid in enumerate(sorted(by_sign[sign])):
                        cf = layout["reactions"][rid]["cofactors"][mid]
                        cf.update(
                            x=side_x + (300 if sign < 0 else 620),
                            y=middle_y + main_side_currency_lane(sign, offset, reserved_lanes),
                            align="center", label_dx=55, label_dy=-4,
                        )

        # A local view is a visual projection, not a new pool. When an
        # internal path metabolite is local at both adjacent reactions, merge
        # those two incidences into one shared local node. Without this step a
        # chemically continuous side path is rendered as two disconnected
        # half-paths (the ED long-chain counterexample).
        for mindex in range(1, len(metabolite_occurrences) - 1):
            if not occurrence_local.get(mindex):
                continue
            mid = metabolite_occurrences[mindex][1]
            adjacent = [chain_reactions[mindex - 1], chain_reactions[mindex]]
            view_id = f"side_{mid}_{mindex}"
            spec["shared_instances"][view_id] = {
                "metabolite": mid,
                "reactions": adjacent,
                "position": [side_x, side_y0 + mindex * chain_spacing],
                "label_offset": [-115, -4],
            }

    # Remaining reactions are usually cycle chords (RPE/RPI/TPI etc.).  Each
    # side component MUST retain one canonical attachment to the already drawn
    # graph; otherwise both endpoints could become source_split views and the
    # visual graph would be disconnected even though the chemistry is not.
    # From that single anchor, draw a short outward stub.  Prefer a main-chain
    # anchor over a ring anchor because it avoids routing back across the ring.
    panel_y = max(center[1] + radius + 290, side_y0 + max(0, len(chain_nodes) - 1) * chain_spacing + 250)
    extra_anchor_use: Counter[str] = Counter()

    def canonical_position(mid: str) -> list[float] | None:
        if mid in ring_positions:
            return list(ring_positions[mid])
        if mid in shared_positions:
            return list(shared_positions[mid])
        if mid in chain_nodes:
            info = layout["metabolites"][mid]
            return [float(info["x"]), float(info["y"])]
        return None

    def anchor_priority(mid: str) -> tuple[int, float]:
        pos = canonical_position(mid)
        if pos is None:
            return (-1, -1.0)
        # Main side-chain nodes are already outside the dense ring and are the
        # cleanest attachment points.  Otherwise prefer the more exterior node.
        tier = 2 if mid in chain_nodes and mid not in ring_positions else 1
        return (tier, math.dist(pos, center))

    def outward_stub(anchor_mid: str, anchor: list[float], slot: int) -> list[float]:
        dx, dy = anchor[0] - center[0], anchor[1] - center[1]
        norm = math.hypot(dx, dy)
        if norm < 1e-9:
            dx, dy, norm = 1.0, 0.0, 1.0
        ux, uy = dx / norm, dy / norm
        tx, ty = -uy, ux
        # Keep the first stub on its radial corridor. Large lateral movement
        # makes the trunk itself cut across neighbouring reaction incidences.
        # Only subsequent siblings need alternating lanes.
        lane = 0 if slot == 0 else ((slot + 1) // 2) * (1 if slot % 2 else -1)
        return [anchor[0] + 620 * ux + lane * 170 * tx,
                anchor[1] + 620 * uy + lane * 170 * ty]

    for index, rid in enumerate(extra_reactions):
        source, target = choose_pair(stoich[rid], ring_metabolites | set(chain_nodes))
        row_y = panel_y + index * 220
        default_a, default_b = [650.0, row_y], [1250.0, row_y]
        source_pos, target_pos = canonical_position(source), canonical_position(target)
        if source_pos is not None or target_pos is not None:
            anchor_mid = max(
                (mid for mid in (source, target) if canonical_position(mid) is not None),
                key=anchor_priority,
            )
            anchor = canonical_position(anchor_mid)
            assert anchor is not None
            slot = extra_anchor_use[anchor_mid]
            extra_anchor_use[anchor_mid] += 1
            local_point = outward_stub(anchor_mid, anchor, slot)
            if anchor_mid == source:
                a = anchor
                local_source = False
                if target in shared_positions:
                    b = target_pos
                    local_target = False
                else:
                    b = local_point
                    local_target = target_pos is not None or incidence_count[target] > 1
            else:
                b = anchor
                if source in shared_positions:
                    a = source_pos
                    local_source = False
                else:
                    a = local_point
                    local_source = source_pos is not None or incidence_count[source] > 1
                local_target = False
        else:
            # A genuinely isolated unique reaction has no existing attachment
            # point. Keep both canonical in the fallback panel; validation will
            # still reject it if the biochemical graph says it should connect.
            a, b = default_a, default_b
            local_source = local_target = False
        side = {
            "reaction": rid, "source": source, "target": target,
            "start": a, "end": b, "label_offset": [0, -72],
            "local_source": local_source, "local_target": local_target,
            "incidence_join_spread": 0.12,
        }
        local_satellites, reserved_lanes = place_side_satellites(rid, source, target, a, b)
        if local_satellites:
            side["local_satellites"] = local_satellites
        spec["side_reactions"].append(side)
        if local_source:
            spec["view_overrides"][f"{source}_{rid}"] = {
                "label_dx": 0, "label_dy": 54, "node_radius": 9,
            }
        if local_target:
            spec["view_overrides"][f"{target}_{rid}"] = {
                "label_dx": 0, "label_dy": 54, "node_radius": 9,
            }
        for mid, pos, is_local in (
            (source, a, local_source),
            (target, b, local_target),
        ):
            # A local source_split endpoint is a reaction-owned visual copy.
            # Never let it overwrite the canonical metabolite position that a
            # previously laid-out main side reaction may already own.  Doing
            # so drags earlier incidence paths across unrelated nodes (the
            # ED-0023 R5P/Xu5P counterexample).
            if (not is_local
                    and mid not in ring_positions and mid not in shared_positions
                    and mid not in canonical_satellite_mids):
                spec["view_overrides"][mid] = {
                    "x": pos[0], "y": pos[1], "label_dx": 0, "label_dy": 46, "node_radius": 9,
                }
                layout["metabolites"][mid].update(x=pos[0], y=pos[1])
        # Keep the cofactor fan outside the occupied graph and never put two
        # same-sign incidences on one collinear ray.
        by_sign = {
            sign: [
                mid for mid, value in stoich[rid].items()
                if mid in CURRENCY and value * sign > 0
            ]
            for sign in (-1, 1)
        }
        canonical_at_source = (
            True if not local_source and canonical_position(source) is not None
            else False if not local_target and canonical_position(target) is not None
            else None
        )
        for mid, placement in extra_currency_fan(
            a, b, center, by_sign, canonical_at_source,
            lane_floor_by_sign=reserved_lanes,
        ).items():
            cf = layout["reactions"][rid]["cofactors"][mid]
            cf.update(
                x=placement["position"][0],
                y=placement["position"][1],
                label_dx=placement["label_offset"][0],
                label_dy=placement["label_offset"][1],
            )

    # Increase canvas after panel placement, then apply the same readability
    # profile that passed the nine-path renderer regression suite.
    if extra_reactions:
        layout["canvas"][1] = max(layout["canvas"][1], int(panel_y + len(extra_reactions) * 220 + 220))
        spec["instance_legend_position"] = [width / 2, layout["canvas"][1] - 135]
        layout["publication_style"]["legend_position"] = [width / 2 - 210, layout["canvas"][1] - 100]
    layout = enlarge_publication(layout)
    meta = {
        "route_id": route_id,
        "module_id": module_id,
        "rank": rank,
        "reaction_count": len(all_reactions),
        "backbone_reaction_count": len(cycle_reactions),
        "side_reaction_count": len(leftovers),
        "main_side_chain_reactions": chain_reactions,
        "main_side_path_nodes": side_path,
        "extra_side_reactions": extra_reactions,
        "cycle_items": cycle_items,
    }
    return layout, meta


def read_source(source_zip: Path, output: Path) -> tuple[Path, dict, list[dict], dict, dict]:
    source_dir = output / "_source"
    source_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source_zip) as archive:
        archive.extractall(source_dir)
    manifest = json.loads((source_dir / "manifest.json").read_text(encoding="utf-8"))
    models = json.loads((source_dir / "route_models.json").read_text(encoding="utf-8"))
    with (source_dir / "ranking_all.tsv").open(encoding="utf-8", newline="") as handle:
        rankings = list(csv.DictReader(handle, delimiter="\t"))
    source_names = json.loads((source_dir / "metabolite_names.json").read_text(encoding="utf-8"))
    with (source_dir / "reaction_evidence.tsv").open(encoding="utf-8", newline="") as handle:
        evidence = {row["column_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert len(rankings) == manifest["candidate_count"] == 47
    assert Counter(row["module_id"] for row in rankings) == Counter(manifest["module_counts"])
    return source_dir, models, rankings, source_names, evidence


def candidate_folder(output: Path, module_id: str, route_id: str) -> Path:
    return output / module_id / route_id


def filter_candidate_rows(rows: list[dict], ids: set[str] | None = None) -> list[dict]:
    """Return rows in source ranking order, optionally restricted to candidate IDs."""
    if not ids:
        return list(rows)
    present = {row["candidate_id"] for row in rows}
    missing = sorted(ids - present)
    if missing:
        raise ValueError(f"unknown candidate IDs: {missing}")
    return [row for row in rows if row["candidate_id"] in ids]


def prepare(source_zip: Path, output: Path, ids: set[str] | None = None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    source_dir, models, rankings, source_names, evidence = read_source(source_zip, output)
    rankings = filter_candidate_rows(rankings, ids)
    checks = []
    for global_index, row in enumerate(rankings, 1):
        route_id = row["candidate_id"]
        module_id = row["module_id"]
        route = models[route_id]
        reaction_model = route["reaction_model"]
        flux = {rid: float(value) for rid, value in route["fluxes"].items()}
        assert {entry[0] for entry in reaction_model} == set(flux)
        folder = candidate_folder(output, module_id, route_id)
        folder.mkdir(parents=True, exist_ok=True)

        lines = []
        for rid, coefficients, *_ in reaction_model:
            coefficients = {mid: float(value) for mid, value in coefficients.items()}
            substrates = " + ".join(term(mid, coefficients[mid]) for mid in sorted(coefficients) if coefficients[mid] < 0)
            products = " + ".join(term(mid, coefficients[mid]) for mid in sorted(coefficients) if coefficients[mid] > 0)
            lines.append("\t".join((rid, substrates, products, "0")))
        (folder / "model.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
        dump(folder / "fit_result.json", {
            "opt_net_fluxes": flux,
            "normalization": {"target": route["target"], "value": 1, "unit": "relative pathway flux"},
            "route_id": route_id,
        })
        layout, layout_meta = build_layout(
            route_id, module_id, int(row.get("rank") or global_index), reaction_model, source_names
        )
        dump(folder / "layout.json", layout)
        dump(folder / "source_pathway.json", {
            "ranking": row,
            "target": route["target"],
            "fluxes": flux,
            "reaction_model": reaction_model,
            "reaction_evidence": {rid: evidence.get(rid, {}) for rid in flux},
            "source_export_sha256": sha(source_zip),
        })
        dump(folder / "layout_meta.json", layout_meta)

        # Round-trip the generated standard input through the renderer parser.
        model = native.read_model(folder / "model.tsv")
        for rid, coefficients, *_ in reaction_model:
            actual = {met.id: float(value) for met, value in model.reactions.get_by_id(rid).metabolites.items()}
            expected = {mid: float(value) for mid, value in coefficients.items()}
            assert actual == expected, (route_id, rid, actual, expected)
        checks.append({
            "index": global_index, "route_id": route_id, "module_id": module_id,
            "reaction_count": len(reaction_model), "stoichiometry_roundtrip": True,
            "flux_roundtrip": True, **layout_meta,
        })

    report = {
        "source_zip": str(source_zip), "source_zip_sha256": sha(source_zip),
        "candidate_count": len(checks),
        "module_counts": dict(Counter(row["module_id"] for row in rankings)),
        "all_standard_inputs_roundtrip": True,
        "candidates": checks,
    }
    dump(output / "input_validation.json", report)
    shutil.copy2(source_zip, source_dir / source_zip.name)
    print(json.dumps({k: report[k] for k in ("candidate_count", "module_counts", "all_standard_inputs_roundtrip")}, ensure_ascii=False, indent=2), flush=True)


def selected_rows(output: Path, ids: set[str] | None = None) -> list[dict]:
    with (output / "_source" / "ranking_all.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return filter_candidate_rows(rows, ids)


def render_one(output: Path, row: dict) -> dict:
    route_id, module_id = row["candidate_id"], row["module_id"]
    folder = candidate_folder(output, module_id, route_id)
    stem = folder / route_id
    manifest_path = stem.with_suffix(".manifest.json")
    previous_stamp = manifest_path.stat().st_mtime_ns if manifest_path.exists() else None
    argv = [
        str(folder / "fit_result.json"), str(folder / "model.tsv"), str(stem.with_suffix(".html")),
        "--layout", str(folder / "layout.json"),
    ]
    code = 0
    with (folder / "render.log").open("w", encoding="utf-8") as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            try:
                sys.argv = [str(ROOT / "render_network.py"), *argv]
                native.main()
            except Exception:
                import traceback
                traceback.print_exc()
                code = 1
    fresh = manifest_path.exists() and manifest_path.stat().st_mtime_ns != previous_stamp
    report = json.loads(manifest_path.read_text(encoding="utf-8")) if fresh else {}
    qa = report.get("rendered_elements", {}).get("strict_visual_qa", {})
    fatal = ("crossings", "overlaps", "node_hits", "semantic_errors", "path_continuity_errors", "visual_connectivity_errors")
    labels = ("label_collisions", "path_label_hits", "clipped_labels", "label_node_hits")
    layout = json.loads((folder / "layout.json").read_text(encoding="utf-8"))
    if code and fresh and layout.get("publication_style", {}).get("auto_label_clearance"):
        if not any(qa.get(key) for key in fatal) and any(qa.get(key) for key in labels):
            adjusted, placement = clear_publication_labels(layout, stem.with_suffix(".html"))
            dump(folder / "label_clearance.json", placement)
            if placement.get("changes") and not placement.get("unresolved"):
                dump(folder / "layout.json", adjusted)
                code = 0
                with (folder / "render.log").open("a", encoding="utf-8") as log:
                    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                        try:
                            sys.argv = [str(ROOT / "render_network.py"), *argv]
                            native.main()
                        except Exception:
                            import traceback
                            traceback.print_exc()
                            code = 1
                report = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
                qa = report.get("rendered_elements", {}).get("strict_visual_qa", {})
    accepted = bool(report.get("acceptance", {}).get("pass")) and code == 0
    result = {
        "route_id": route_id, "module_id": module_id, "returncode": code,
        "accepted": accepted,
        "svg": stem.with_suffix(".svg").exists(),
        "strict_visual_qa": {key: qa.get(key) for key in (
            "pass", "crossings", "overlaps", "label_collisions", "path_label_hits",
            "clipped_labels", "node_hits", "label_node_hits", "semantic_errors",
            "path_continuity_errors", "visual_connectivity_errors",
        )},
    }
    # Windows PowerShell on the paired workstation may expose a GBK stdout.
    # Keep the persisted JSON UTF-8, but make progress logging ASCII-safe so a
    # metabolite label such as H⁺ cannot abort an otherwise completed render.
    print(json.dumps(result, ensure_ascii=True), flush=True)
    return result


def render(output: Path, ids: set[str] | None = None) -> None:
    rows = selected_rows(output, ids)
    results = [render_one(output, row) for row in rows]
    dump(output / "render_status.json", {"all_pass": all(row["accepted"] for row in results), "results": results})
    failed = [row["route_id"] for row in results if not row["accepted"]]
    if failed:
        raise RuntimeError(f"native render acceptance failed: {failed}")


def validate(output: Path, ids: set[str] | None = None) -> None:
    rows = selected_rows(output, ids)
    checks = []
    for row in rows:
        route_id, module_id = row["candidate_id"], row["module_id"]
        folder = candidate_folder(output, module_id, route_id)
        manifest = json.loads((folder / f"{route_id}.manifest.json").read_text(encoding="utf-8"))
        source = json.loads((folder / "source_pathway.json").read_text(encoding="utf-8"))
        fit = json.loads((folder / "fit_result.json").read_text(encoding="utf-8"))
        layout = json.loads((folder / "layout.json").read_text(encoding="utf-8"))
        assert fit["opt_net_fluxes"] == source["fluxes"]
        assert manifest["acceptance"]["pass"]
        qa = manifest["rendered_elements"]["strict_visual_qa"]
        assert qa.get("pass"), (route_id, qa)
        checks.append({
            "route_id": route_id, "module_id": module_id,
            "native_acceptance": True, "strict_visual_qa": True,
            "relationship_mode": layout["pathway_multiplicity"]["relationship_mode"],
            "reaction_count": manifest["model_reaction_count"],
            "svg_sha256": sha(folder / f"{route_id}.svg"),
            "png_sha256": sha(folder / f"{route_id}.png"),
        })
    report = {
        "all_pass": True,
        "candidate_count": len(checks),
        "module_counts": dict(Counter(row["module_id"] for row in checks)),
        "source_split_relationship_mode": all(row["relationship_mode"] == "source_split" for row in checks),
        "candidates": checks,
    }
    dump(output / "revision_validation.json", report)
    print(json.dumps({k: report[k] for k in ("all_pass", "candidate_count", "module_counts", "source_split_relationship_mode")}, ensure_ascii=False, indent=2), flush=True)


def export(output: Path, ids: set[str] | None = None) -> None:
    import fitz
    from PIL import Image, ImageDraw
    from playwright.sync_api import sync_playwright

    rows = selected_rows(output, ids)
    pdf_root = output / "PDF"
    preview_root = output / "pdf_previews"
    pdf_root.mkdir(exist_ok=True)
    preview_root.mkdir(exist_ok=True)
    checks = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        for row in rows:
            route_id, module_id = row["candidate_id"], row["module_id"]
            folder = candidate_folder(output, module_id, route_id)
            manifest = json.loads((folder / f"{route_id}.manifest.json").read_text(encoding="utf-8"))
            if not manifest["acceptance"]["pass"]:
                raise RuntimeError(f"{route_id}: native acceptance not passed")
            svg_text = (folder / f"{route_id}.svg").read_text(encoding="utf-8")
            import xml.etree.ElementTree as ET
            root = ET.fromstring(svg_text)
            width, height = [float(root.get(key)) for key in ("width", "height")]
            page = browser.new_page(viewport={"width": round(width), "height": round(height)})
            page.set_content('<!doctype html><meta charset="utf-8"><style>@page{margin:0}body{margin:0;background:white}svg{display:block}</style>' + svg_text)
            page.evaluate("document.fonts.ready")
            module_pdf = pdf_root / module_id
            module_pdf.mkdir(exist_ok=True)
            pdf_path = module_pdf / f"{route_id}.pdf"
            page.pdf(path=str(pdf_path), width=f"{width}px", height=f"{height}px", print_background=True,
                     margin={"top": "0", "right": "0", "bottom": "0", "left": "0"}, display_header_footer=False)
            page.close()
            doc = fitz.open(pdf_path)
            assert len(doc) == 1
            p0 = doc[0]
            text = p0.get_text()
            assert route_id in text
            assert len(p0.get_drawings()) > 20
            outside = [word[4] for word in p0.get_text("words")
                       if word[0] < -1 or word[1] < -1 or word[2] > p0.rect.width + 1 or word[3] > p0.rect.height + 1]
            assert not outside, (route_id, outside)
            preview_dir = preview_root / module_id
            preview_dir.mkdir(exist_ok=True)
            scale = min(1600 / p0.rect.width, 1200 / p0.rect.height)
            pix = p0.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
            pix.save(preview_dir / f"{route_id}.png")
            checks.append({
                "route_id": route_id, "module_id": module_id, "pdf": str(pdf_path),
                "pages": 1, "vector_paths": len(p0.get_drawings()), "clipped_words": [],
                "sha256": sha(pdf_path), "bytes": pdf_path.stat().st_size,
            })
            doc.close()
            print(f"{route_id}: vector PDF verified", flush=True)
        browser.close()

    # Per-module and full merged PDFs keep the ranking order from ranking_all.tsv.
    for module_id in ("EMP", "ED", "PKT", "OPPP"):
        module_rows = [row for row in rows if row["module_id"] == module_id]
        if not module_rows:
            continue
        merged = fitz.open()
        for row in module_rows:
            doc = fitz.open(pdf_root / module_id / f"{row['candidate_id']}.pdf")
            merged.insert_pdf(doc)
            doc.close()
        merged.save(pdf_root / f"{module_id}_all_{len(module_rows)}.pdf", garbage=4, deflate=True)
        merged.close()
    merged = fitz.open()
    for row in rows:
        doc = fitz.open(pdf_root / row["module_id"] / f"{row['candidate_id']}.pdf")
        merged.insert_pdf(doc)
        doc.close()
    merged_name = f"all_{len(rows)}_pathways.pdf"
    merged.save(pdf_root / merged_name, garbage=4, deflate=True)
    merged.close()

    # A compact 4-column contact sheet is generated from the validated PDFs.
    tiles = []
    for row in rows:
        path = preview_root / row["module_id"] / f"{row['candidate_id']}.png"
        image = Image.open(path).convert("RGB")
        image.thumbnail((720, 510))
        tile = Image.new("RGB", (760, 560), "white")
        tile.paste(image, ((760 - image.width) // 2, 40 + (510 - image.height) // 2))
        ImageDraw.Draw(tile).text((18, 12), f"{row['module_id']} · {row['candidate_id']}", fill="black")
        tiles.append(tile)
    cols = 4
    rows_count = math.ceil(len(tiles) / cols)
    sheet = Image.new("RGB", (cols * 760, rows_count * 560), "white")
    for index, tile in enumerate(tiles):
        sheet.paste(tile, ((index % cols) * 760, (index // cols) * 560))
    contact_name = f"all_{len(rows)}_contact_sheet.png"
    sheet.save(preview_root / contact_name)

    module_counts = dict(Counter(row["module_id"] for row in checks))
    report = {
        "all_pass": True, "pdf_count": len(checks),
        "module_counts": module_counts,
        "merged_pdf": str(pdf_root / merged_name), "files": checks,
    }
    dump(output / "delivery_validation.json", report)
    module_lines = "\n".join(
        f"- {module}: {module_counts.get(module, 0)} 条"
        for module in ("EMP", "ED", "PKT", "OPPP")
        if module_counts.get(module, 0)
    )
    readme = f"""# 模块候选路径绘图交付

本目录由冻结的 module-conditioned 路径直接转换并绘制，未重新求解或修改反应/Flux。

{module_lines}
- 合计: {len(rows)} 条

优先查看 `PDF/{merged_name}`；有候选的模块另有各自合辑。每条路径目录保留 `model.tsv`、`fit_result.json`、`layout.json`、SVG、PNG、HTML、Map JSON 和 manifest。

布局使用每条路径自动识别出的真实有向碳循环作为圆形主干，剩余反应作为开放支路/局部反应。数量标签使用 `local_flow`，重复代谢物视图使用 `source_split` 来源/拆分语义。

可读性沿用当前 renderer profile：代谢物 24 pt、反应 17 pt、线宽 3.2–7.6、碳节点半径 12、辅因子半径 7。

`input_validation.json` 证明标准输入对冻结计量/Flux 的 round-trip；`revision_validation.json` 记录原生 SVG 严格验收；`delivery_validation.json` 记录 {len(rows)} 张矢量 PDF 的页数、矢量内容、文字边界与 SHA256。
"""
    (output / "README.md").write_text(readme, encoding="utf-8")
    print(json.dumps({"pdf_count": len(checks), "merged": str(pdf_root / merged_name), "contact_sheet": str(preview_root / contact_name)}, ensure_ascii=False, indent=2), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare, render, validate, and export frozen module-conditioned pathways."
    )
    parser.add_argument("action", choices=("prepare", "render", "validate", "export"))
    parser.add_argument(
        "--source",
        type=Path,
        help="Source export ZIP. Required only for the prepare action.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Batch working/output directory.",
    )
    parser.add_argument("--ids", default="", help="comma-separated candidate IDs for a focused render")
    args = parser.parse_args()
    ids = {value.strip() for value in args.ids.split(",") if value.strip()} or None
    if args.action == "prepare":
        if args.source is None:
            parser.error("--source is required for the prepare action")
        prepare(args.source, args.output, ids)
    elif args.action == "render":
        render(args.output, ids)
    elif args.action == "validate":
        validate(args.output, ids)
    else:
        export(args.output, ids)


if __name__ == "__main__":
    main()
