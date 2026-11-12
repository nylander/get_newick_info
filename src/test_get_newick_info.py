import math
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PROGRAM = ROOT / "src" / "get_newick_info"
TAXA = (
    "Andricus_curvator\n"
    "Andricus_kollari\n"
    "Barbotinia_oraniensis\n"
    "Diplolepis_rosae\n"
    "Panteliella_bicolor\n"
    "Periclistus_brandtii\n"
    "Plagiotrochus_quercusilicis\n"
    "Synergus_gallaepomiformis\n"
)
SUPPORT = (
    "0.999\n1\n0.60099999999999998\n"
    "0.58899999999999997\n0.47699999999999998\n"
)
TL = "0.48612273009999996\n"


class OutputTests(unittest.TestCase):
    def run_info(self, *args):
        result = subprocess.run(
            [str(PROGRAM), *map(str, args)],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(result.stderr, "")
        return result.stdout

    def summary_metric(self, output, label):
        prefix = f"{label}: "
        return float(next(line[len(prefix):] for line in output.splitlines()
                          if line.startswith(prefix)))

    def test_value_only_options(self):
        for option, expected in (
            ("--taxa", TAXA),
            ("-t", TAXA),
            ("--support", SUPPORT),
            ("--tl", "0.48455999999999999\n"),
        ):
            with self.subTest(option=option):
                self.assertEqual(
                    self.run_info(option, ROOT / "data" / "brlens_pp.tree"),
                    expected,
                )

    def test_output_options_are_mutually_exclusive(self):
        input_path = ROOT / "data" / "brlens_pp.tree"
        for options in (
            ("--taxa", "--support"),
            ("--taxa", "--tl"),
            ("--support", "--tl"),
            ("-t", "--support"),
        ):
            with self.subTest(options=options):
                result = subprocess.run(
                    [str(PROGRAM), *options, str(input_path)],
                    capture_output=True,
                    text=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, "")
                self.assertIn("mutually exclusive", result.stderr)

    def test_tree_length_example(self):
        self.assertEqual(
            self.run_info("--tl", ROOT / "data" / "brlens.tree"), TL
        )

    def test_multi_tree_and_multi_file_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            path.write_text("(B:1,A:2)90;(A:3,C:4)80;\n")
            self.assertEqual(self.run_info("--tl", path, path), "3\n7\n3\n7\n")
            self.assertEqual(self.run_info("--taxa", path), "A\nB\nC\n")
            self.assertEqual(self.run_info("--support", path), "90\n80\n")
            self.assertIn("Tree length (TL): 10\n", self.run_info(path))
            path.write_text("(A,B);\n")
            self.assertEqual(self.run_info("--tl", path), "0\n")
            self.assertEqual(self.run_info("--support", path), "")
            path.write_text("")
            self.assertEqual(self.run_info("--tl", path), "")
            path.write_text("(A:1,B:2);" * 17)
            self.assertEqual(self.run_info("--tl", path), "3\n" * 17)

    def test_tl_does_not_print_partial_results_on_parse_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.nwk"
            path.write_text("(A:1,B:2);(C:3,D:4)")
            result = subprocess.run(
                [str(PROGRAM), "--tl", str(path)], capture_output=True, text=True
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(result.stdout, "")
            self.assertIn("expected ';'", result.stderr)

    def test_default_summary(self):
        path = ROOT / "data" / "brlens.tree"
        output = self.run_info(path)
        self.assertEqual(
            output,
            f"File: {path}\n"
            "Number of trees: 1\n"
            "Number of unique tip labels: 8\n"
            "Branch lengths present: true\n"
            "Rooted tree: false\n"
            "Dichotomous: false\n"
            "Comments present (NHX or Enhanced Newick format): false\n"
            "Is the tree ultrametric: false\n"
            "Tree height: NA\n"
            f"Tree length (TL): {TL}"
            "Mean resampling proportion (MRP): 0.54600000000000004\n"
            "Standard deviation of resampling proportion (s_RP): 0.25726251184344756\n"
            "Mean character support (MCS): 1.4038482366507115\n"
            "Standard deviation of character support (s_CS): 1.6167469061723576\n"
            "Sum of support values: 274\n\n",
        )

    def test_mrp_and_mcs_from_resampling_proportions(self):
        output = self.run_info(ROOT / "data" / "brlens_pp.tree")
        supports = [0.999, 0.99, 0.601, 0.589, 0.477]
        expected_mrp = sum(supports) / 5
        expected_mcs = sum(-math.log1p(-support) for support in supports) / 5
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean resampling proportion (MRP)"),
            expected_mrp,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean character support (MCS)"),
            expected_mcs,
        )
        cs = [-math.log1p(-support) for support in supports]
        rp_mean = sum(supports) / len(supports)
        cs_mean = sum(cs) / len(cs)

        self.assertAlmostEqual(
            self.summary_metric(
                output, "Standard deviation of resampling proportion (s_RP)"
            ),
            math.sqrt(sum((x - rp_mean) ** 2 for x in supports) / len(supports)),
        )
        self.assertAlmostEqual(
            self.summary_metric(
                output, "Standard deviation of character support (s_CS)"
            ),
            math.sqrt(sum((x - cs_mean) ** 2 for x in cs) / len(cs)),
        )

    def test_supports_on_zero_to_100_scale(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "supports.nwk"
            path.write_text("((A,B)100,(C,D)99,E);\n")
            output = self.run_info(path)
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean resampling proportion (MRP)"),
            0.99,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean character support (MCS)"),
            -math.log(0.01),
        )

    def test_support_metrics_use_available_supports(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "partial-supports.nwk"
            path.write_text("((A,B)0.5,C,D,E);\n")
            output = self.run_info(path)
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean resampling proportion (MRP)"),
            0.25,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean character support (MCS)"),
            -math.log(0.5) / 2,
        )
        self.assertEqual(
            self.summary_metric(output, "Standard deviation of resampling proportion (s_RP)"),
            0.0,
        )
        self.assertEqual(
            self.summary_metric(output, "Standard deviation of character support (s_CS)"),
            0.0,
        )

    def test_support_metrics_are_na_without_supports(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unsupported.nwk"
            path.write_text("(A,B,C,D);\n")
            output = self.run_info(path)
        self.assertIn("Mean resampling proportion (MRP): NA\n", output)
        self.assertIn("Mean character support (MCS): NA\n", output)
        self.assertIn("Standard deviation of resampling proportion (s_RP): NA\n", output)
        self.assertIn("Standard deviation of character support (s_CS): NA\n", output)

    def test_support_metrics_are_na_when_tip_count_has_no_denominator(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small-supported.nwk"
            path.write_text("((A,B)0.5,C);\n")
            output = self.run_info(path)
        self.assertIn("Mean resampling proportion (MRP): NA\n", output)
        self.assertIn("Mean character support (MCS): NA\n", output)
        self.assertIn("Standard deviation of resampling proportion (s_RP): NA\n", output)
        self.assertIn("Standard deviation of character support (s_CS): NA\n", output)

    def test_multiple_tree_metrics_are_averaged(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            path.write_text("((A,B)0.5,C,D);((E,F)1.0,G,H);\n")
            output = self.run_info(path)
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean resampling proportion (MRP)"),
            0.745,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "Mean character support (MCS)"),
            (-math.log(0.5) - math.log(0.01)) / 2,
        )

    def test_multiple_tree_standard_deviations_are_averaged(self):
        trees = (
            ("((A,B)0.2,(C,D)0.8,E);", [0.2, 0.8]),
            ("((A,B)20,(C,D)40,(E,F)60,G);", [0.2, 0.4, 0.6]),
            ("((A,B)1.0,C,D);", [0.99]),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            path.write_text("".join(tree for tree, _ in trees))
            output = self.run_info(path)

        for label, transform in (
            ("Standard deviation of resampling proportion (s_RP)", lambda rp: rp),
            ("Standard deviation of character support (s_CS)", lambda rp: -math.log1p(-rp)),
        ):
            deviations = []
            for _, supports in trees:
                values = [transform(rp) for rp in supports]
                mean = sum(values) / len(values)
                deviations.append(
                    math.sqrt(sum((value - mean) ** 2 for value in values) / len(values))
                )
            with self.subTest(label=label):
                self.assertAlmostEqual(
                    self.summary_metric(output, label),
                    sum(deviations) / len(deviations),
                )

    def test_multiple_tree_standard_deviations_are_na_if_any_tree_is_unsupported(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            for trees in (
                "((A,B)0.2,(C,D)0.8,E);(A,B,C,D);",
                "(A,B,C,D);((A,B)0.2,(C,D)0.8,E);",
            ):
                with self.subTest(trees=trees):
                    path.write_text(trees)
                    output = self.run_info(path)
                    self.assertIn(
                        "Standard deviation of resampling proportion (s_RP): NA\n", output
                    )
                    self.assertIn(
                        "Standard deviation of character support (s_CS): NA\n", output
                    )


if __name__ == "__main__":
    unittest.main()
