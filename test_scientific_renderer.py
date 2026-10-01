import copy
import json
import math
import unittest
from pathlib import Path

from scripts_render_d3flux import (
    apply_initial_map,
    apply_layout,
    apply_line_style,
    read_model,
)
from scientific_layout import apply_scientific_layout, build_publication_qa


ROOT = Path(__file__).resolve().parent
DEMO = ROOT / "candidate01_demo"


class ScientificRendererContractTests(unittest.TestCase):
    def setUp(self):
        self.layout = json.loads((DEMO / "layout.json").read_text(encoding="utf-8"))
        # Retain regression coverage of the existing legacy layout mode.
        # The default circular mode has its own stricter contract suite.
        self.layout.pop('circular_backbone', None)
        self.layout['canvas'] = [1600, 1000]
        self.fit = json.loads((DEMO / "fit_result.json").read_text(encoding="utf-8"))

    def fresh_model(self):
        return read_model(DEMO / "model.tsv")

    def test_p0_default_initial_map_cannot_override_coordinates(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        expected = {
            obj.id: (obj.notes["map_info"]["x"], obj.notes["map_info"]["y"])
            for obj in list(model.metabolites) + list(model.reactions)
            if not obj.notes["map_info"].get("hidden")
        }
        apply_initial_map(model, ROOT / "initial_template.json", inherit_coordinates=False)
        actual = {
            obj.id: (obj.notes["map_info"]["x"], obj.notes["map_info"]["y"])
            for obj in list(model.metabolites) + list(model.reactions)
            if not obj.notes["map_info"].get("hidden")
        }
        self.assertEqual(expected, actual)

    def test_p0_coordinate_inheritance_requires_explicit_opt_in(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        self.assertEqual((493, 170), (
            model.metabolites.get_by_id("F6P").notes["map_info"]["x"],
            model.metabolites.get_by_id("F6P").notes["map_info"]["y"],
        ))
        apply_initial_map(model, ROOT / "initial_template.json", inherit_coordinates=True)
        inherited = model.metabolites.get_by_id("F6P").notes["map_info"]
        self.assertNotEqual((493, 170), (inherited["x"], inherited["y"]))

    def test_p1_scientific_layout_assigns_roles_and_cycle_skeleton(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        summary = apply_scientific_layout(model, self.layout)
        self.assertTrue(summary["enabled"])
        self.assertEqual("source", model.metabolites.get_by_id("Formaldehyde").notes["map_info"]["scientific_role"])
        self.assertEqual("sink", model.metabolites.get_by_id("Pyruvate").notes["map_info"]["scientific_role"])
        self.assertEqual("cycle_anchor", model.metabolites.get_by_id("Ru5P").notes["map_info"]["scientific_role"])
        self.assertEqual(0, summary["lanes"])
        structures = {item["name"]: item for item in summary["structures"]}
        self.assertIn("carbon_cycle_arc", structures)
        self.assertIn("ppp_regeneration_ring", structures)
        outer = structures["carbon_cycle_arc"]["items"]
        self.assertEqual("Ru5P", outer[0])
        self.assertEqual("G3P", outer[-1])
        ys = []
        for item_id in outer:
            collection = model.reactions if item_id in model.reactions else model.metabolites
            info = collection.get_by_id(item_id).notes["map_info"]
            self.assertEqual("carbon_cycle_arc", info["layout_structure"])
            ys.append(round(float(info["y"]), 3))
        self.assertGreater(len(set(ys)), 5)  # regression: backbone must not collapse to one horizontal lane

        order = summary["structure_order_optimization"]
        self.assertEqual(1, len(order))
        self.assertTrue(order[0]["optimized"])
        self.assertNotEqual(order[0]["order_before"], order[0]["order_after"])

    def test_p2_cofactors_are_local_satellites_not_global_layout_nodes(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        for mid in ("NADP", "NADPH", "H", "H2O"):
            self.assertTrue(model.metabolites.get_by_id(mid).notes["map_info"].get("hidden"))
        g6 = model.reactions.get_by_id("IML_G6PDH2r").notes["map_info"]["cofactors"]
        self.assertEqual({"NADP", "NADPH", "H"}, set(g6))

    def test_p3_publication_mode_uses_width_not_low_flux_dashes(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        flux = {rid: float(value) for rid, value in self.fit["opt_net_fluxes"].items()}
        result = apply_line_style(model, flux, self.layout["line_style"])
        self.assertEqual("flux_width_only", result["mode"])
        self.assertEqual([], result["low_relative_flux_reactions"])
        widths = {r.id: r.notes["map_info"]["stroke"] for r in model.reactions}
        self.assertGreater(widths["RUMP_R05338"], widths["IML_PGI"])
        self.assertTrue(all(r.notes["map_info"]["dasharray"] == "none" for r in model.reactions))

    def test_publication_qa_passes_candidate01_scene(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        apply_scientific_layout(model, self.layout)
        qa = build_publication_qa(model, self.layout)
        self.assertTrue(qa["pass"], qa)
        self.assertFalse(qa["missing_coordinates"])
        self.assertFalse(qa["out_of_bounds"])
        self.assertTrue(all(x["strictly_monotonic"] for x in qa["lane_checks"]))
        self.assertTrue(all(x["pass"] for x in qa["structure_checks"]), qa)
        self.assertEqual(0, qa["route_visual"]["near_edge_count"], qa)
        self.assertEqual(0, qa["route_visual"]["route_stretch_violation_count"], qa)

    def test_p4_route_optimizer_does_not_leave_candidate_crossings(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        summary = apply_scientific_layout(model, self.layout)
        routes = summary["routes"]
        self.assertEqual(0, routes["estimated_crossings"], routes.get("crossing_pairs"))
        self.assertEqual(0, routes["estimated_obstacle_hits"], routes.get("obstacle_hits"))
        bundle = routes["bundles"]["ppp_regeneration"]["optimization"]
        self.assertTrue(bundle["optimized"])
        self.assertLessEqual(bundle["crossings_after"], bundle["crossings_before"])
        edge = routes["edge_optimization"]
        self.assertTrue(edge["enabled"])
        self.assertLessEqual(edge["crossings_after"], edge["crossings_before"])
        hub = routes["hub_optimization"]
        self.assertTrue(hub["enabled"])
        self.assertLessEqual(hub["crossings_after"], hub["crossings_before"])
        waypoint = routes["waypoint_optimization"]
        self.assertTrue(waypoint["enabled"])
        self.assertEqual(0, waypoint["crossings_after"])
        node = routes["node_optimization"]
        self.assertTrue(node["enabled"])
        self.assertEqual(0, node["crossings_after"])
        self.assertEqual(0, node["obstacles_after"])
        compaction = routes["compaction"]
        self.assertTrue(compaction["enabled"])
        self.assertEqual(0, routes["estimated_crossings"])
        self.assertEqual(0, routes["estimated_obstacle_hits"])
        if compaction["applied"]:
            self.assertLessEqual(compaction["route_length_after"],
                                 compaction["route_length_before"])

    def test_p5_multi_endpoint_reactions_are_grouped_hyperedges(self):
        model = self.fresh_model()
        apply_layout(model, self.layout)
        summary = apply_scientific_layout(model, self.layout)
        expected = {
            "RUMP_R05338",
            "IML_EDA",
            "IML_TKT2",
            "IML_TALA",
            "IML_TKT1",
        }
        self.assertEqual(len(expected), summary["routes"]["reaction_buses"])
        for reaction_id in expected:
            reaction = model.reactions.get_by_id(reaction_id)
            info = reaction.notes["map_info"]
            self.assertEqual("reaction_bus", info["route_style"])
            bus = info["reaction_bus"]
            self.assertGreater(float(bus["half_length"]), 0)
            self.assertEqual(2, len(bus["axis"]))
            self.assertTrue(info["edge_routes"])
            for edge in info["edge_routes"].values():
                self.assertEqual("reaction_bus", edge["style"])
                self.assertEqual(0.0, edge["channel_offset"])
                self.assertNotIn("waypoint_before", edge)
                self.assertNotIn("waypoint_after", edge)

        # Single-substrate/single-product reactions keep ordinary smooth
        # routing; the hyperedge rule must not flatten the whole map.
        self.assertNotEqual(
            "reaction_bus",
            model.reactions.get_by_id("IML_PGI").notes["map_info"]["route_style"],
        )


if __name__ == "__main__":
    unittest.main()
