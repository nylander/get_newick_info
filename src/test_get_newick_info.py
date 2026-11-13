import math
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PROGRAM = ROOT / "src" / "get_newick_info"
EXAMPLE = ROOT / "data" / "brlens_bp.tree"
HEADER = (
    "#\tNtips\tblens\tsupport\troot\tdicho\tcomm\tultra\t"
    "height\tTL\tsumSup\tMRP\tsRP\tMCS\tsCS\tNtrees\tfile\n"
)
EXAMPLE_ROW = (
    "1\t8\ttrue\ttrue\tfalse\tfalse\tfalse\tfalse\tNA\t"
    "0.48529400210000007\t288\t0.57400000000000007\t"
    "0.24368832553078942\t1.4543634882709771\t1.5904965740184576\t"
    "1\tbrlens_bp.tree\n"
)
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
PP_SUPPORT = (
    "0.999\n1\n0.60099999999999998\n"
    "0.58899999999999997\n0.47699999999999998\n"
)
BRANCH_LENGTHS = (
    "0.0078347836000000008\n0.028452826099999998\n"
    "0.025623778400000002\n0.058971544100000002\n0.0069431998\n"
    "0.0644708066\n0.0036897742999999999\n0.1440215645\n"
    "0.024767003999999999\n0.017712700099999999\n"
    "0.0086275608999999993\n0.062971422999999999\n0.031207036699999999\n"
)


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

    def summary(self, output):
        summaries = self.summaries(output)
        self.assertEqual(len(summaries), 1)
        return summaries[0]

    def summaries(self, output):
        lines = output.splitlines()
        self.assertTrue(lines)
        self.assertEqual(lines[0] + "\n", HEADER)
        columns = lines[0].split("\t")
        reports = []
        for line in lines[1:]:
            values = line.split("\t")
            self.assertEqual(len(values), len(columns))
            reports.append(dict(zip(columns, values)))
        return reports

    def summary_metric(self, output, label):
        return float(self.summary(output)[label])

    def test_value_only_options(self):
        for option, expected in (
            ("--taxa", TAXA),
            ("-t", TAXA),
            ("--support", "60\n23\n52\n53\n100\n"),
            ("-s", "60\n23\n52\n53\n100\n"),
            ("--brlens", BRANCH_LENGTHS),
            ("-b", BRANCH_LENGTHS),
            ("--tl", "0.48529400210000007\n"),
        ):
            with self.subTest(option=option):
                self.assertEqual(self.run_info(option, EXAMPLE), expected)
        self.assertEqual(
            self.run_info("--support", ROOT / "data" / "brlens_pp.tree"),
            PP_SUPPORT,
        )

    def test_output_options_are_mutually_exclusive(self):
        input_path = ROOT / "data" / "brlens_pp.tree"
        for options in (
            ("--taxa", "--support"),
            ("--taxa", "--tl"),
            ("--support", "--tl"),
            ("-t", "--support"),
            ("--brlens", "--taxa"),
            ("--brlens", "--support"),
            ("--brlens", "--tl"),
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
            self.run_info("--tl", ROOT / "data" / "brlens.tree"),
            "0.48455999999999999\n",
        )

    def test_multi_tree_and_multi_file_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            path.write_text("(B:1,A:2)90;(A:3,C:4)80;\n")
            self.assertEqual(self.run_info("--tl", path, path), "3\n7\n3\n7\n")
            self.assertEqual(self.run_info("--taxa", path), "A\nB\nC\n")
            self.assertEqual(self.run_info("--support", path), "90\n80\n")
            reports = self.summaries(self.run_info(path))
            self.assertEqual([report["TL"] for report in reports], ["3", "7"])
            self.assertEqual([report["#"] for report in reports], ["1", "2"])
            self.assertEqual([report["Ntrees"] for report in reports], ["2", "2"])
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

    def test_branch_lengths_across_trees_and_files(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            for option in ("--brlens", "-b"):
                with self.subTest(option=option):
                    path.write_text("((A:1,B:2)90:3,C:4):5;(D:6,E:7);")
                    expected = "1\n2\n3\n4\n5\n6\n7\n"
                    self.assertEqual(self.run_info(option, path, path), expected * 2)
                    path.write_text("(A:1,B:2);" * 17)
                    self.assertEqual(self.run_info(option, path), "1\n2\n" * 17)
                    path.write_text("(A,B);")
                    self.assertEqual(self.run_info(option, path), "")
                    path.write_text("")
                    self.assertEqual(self.run_info(option, path), "")

    def test_branch_lengths_do_not_print_partial_results_on_parse_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.nwk"
            for trees, error in (
                ("(A:1,B:2);(C:3,D:4)", "expected ';'"),
                ("(A:1,B:2);(C:3,)", "expected a taxon label"),
                ("(C:3,)", "expected a taxon label"),
            ):
                with self.subTest(trees=trees):
                    path.write_text(trees)
                    result = subprocess.run(
                        [str(PROGRAM), "--brlens", str(path)],
                        capture_output=True,
                        text=True,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertEqual(result.stdout, "")
                    self.assertIn(error, result.stderr)

    def test_default_summary(self):
        self.assertEqual(self.run_info(EXAMPLE), HEADER + EXAMPLE_ROW)

    def test_multi_tree_file_reports_each_tree(self):
        expected = (
            "1\t8\ttrue\ttrue\tfalse\tfalse\tfalse\tfalse\tNA\t"
            "0.48529400210000007\t288\t0.57400000000000007\t"
            "0.24368832553078942\t1.4543634882709771\t1.5904965740184576\t"
            "3\tmulti.trees\n"
            "2\t8\ttrue\ttrue\tfalse\tfalse\tfalse\tfalse\tNA\t"
            "0.4637315622\t298\t0.59599999999999997\t"
            "0.2218648237102944\t1.1016457698187057\t0.6829615010376946\t"
            "3\tmulti.trees\n"
            "3\t8\ttrue\ttrue\tfalse\tfalse\tfalse\tfalse\tNA\t"
            "0.47744841850000003\t308\t0.6140000000000001\t"
            "0.24630062931304095\t1.5586228729976344\t1.5546836427319075\t"
            "3\tmulti.trees\n"
        )
        self.assertEqual(
            self.run_info(ROOT / "data" / "multi.trees"),
            HEADER + expected,
        )

    def test_summary_output_options(self):
        for option in ("-n", "--noverbose"):
            with self.subTest(option=option):
                self.assertEqual(self.run_info(option, EXAMPLE), EXAMPLE_ROW)
        fullpath_row = EXAMPLE_ROW.replace("brlens_bp.tree\n", f"{EXAMPLE.resolve()}\n")
        for option in ("-p", "--fullpath"):
            with self.subTest(option=option):
                self.assertEqual(self.run_info(option, EXAMPLE), HEADER + fullpath_row)
        self.assertEqual(self.run_info("-n", "-p", EXAMPLE), fullpath_row)

    def test_multiple_file_summaries(self):
        self.assertEqual(
            self.run_info(EXAMPLE, EXAMPLE),
            (HEADER + EXAMPLE_ROW) * 2,
        )
        self.assertEqual(self.run_info("-n", EXAMPLE, EXAMPLE), EXAMPLE_ROW * 2)

    def test_summary_numeric_options(self):
        output = self.summary(self.run_info("--lp", "0.1", EXAMPLE))
        self.assertAlmostEqual(float(output["MRP"]), (0.6 + 0.23 + 0.52 + 0.53 + 0.9) / 5)
        self.assertAlmostEqual(
            float(output["MCS"]),
            sum(-math.log1p(-rp) for rp in (0.6, 0.23, 0.52, 0.53, 0.9)) / 5,
        )
        path = ROOT / "data" / "ultrametric.tree"
        loose = self.summary(self.run_info("--tol", "1e-7", path))
        strict = self.summary(self.run_info("--tol", "1e-11", path))
        self.assertEqual((loose["ultra"], loose["height"]), ("true", "1"))
        self.assertEqual((strict["ultra"], strict["height"]), ("false", "NA"))

    def test_mrp_and_mcs_from_resampling_proportions(self):
        output = self.run_info(ROOT / "data" / "brlens_pp.tree")
        supports = [0.999, 0.99, 0.601, 0.589, 0.477]
        expected_mrp = sum(supports) / 5
        expected_mcs = sum(-math.log1p(-support) for support in supports) / 5
        self.assertAlmostEqual(
            self.summary_metric(output, "MRP"),
            expected_mrp,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "MCS"),
            expected_mcs,
        )
        cs = [-math.log1p(-support) for support in supports]
        rp_mean = sum(supports) / len(supports)
        cs_mean = sum(cs) / len(cs)

        self.assertAlmostEqual(
            self.summary_metric(
                output, "sRP"
            ),
            math.sqrt(sum((x - rp_mean) ** 2 for x in supports) / len(supports)),
        )
        self.assertAlmostEqual(
            self.summary_metric(
                output, "sCS"
            ),
            math.sqrt(sum((x - cs_mean) ** 2 for x in cs) / len(cs)),
        )

    def test_supports_on_zero_to_100_scale(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "supports.nwk"
            path.write_text("((A,B)100,(C,D)99,E);\n")
            output = self.run_info(path)
        self.assertAlmostEqual(
            self.summary_metric(output, "MRP"),
            0.99,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "MCS"),
            -math.log(0.01),
        )

    def test_support_metrics_use_available_supports(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "partial-supports.nwk"
            path.write_text("((A,B)0.5,C,D,E);\n")
            output = self.run_info(path)
        self.assertAlmostEqual(
            self.summary_metric(output, "MRP"),
            0.25,
        )
        self.assertAlmostEqual(
            self.summary_metric(output, "MCS"),
            -math.log(0.5) / 2,
        )
        self.assertEqual(
            self.summary_metric(output, "sRP"),
            0.0,
        )
        self.assertEqual(
            self.summary_metric(output, "sCS"),
            0.0,
        )

    def test_support_metrics_are_na_without_supports(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unsupported.nwk"
            path.write_text("(A,B,C,D);\n")
            output = self.run_info(path)
        report = self.summary(output)
        for column in ("MRP", "sRP", "MCS", "sCS"):
            self.assertEqual(report[column], "NA")

    def test_support_metrics_are_na_when_tip_count_has_no_denominator(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small-supported.nwk"
            path.write_text("((A,B)0.5,C);\n")
            output = self.run_info(path)
        report = self.summary(output)
        for column in ("MRP", "sRP", "MCS", "sCS"):
            self.assertEqual(report[column], "NA")

    def test_multiple_tree_metrics_are_reported_per_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            path.write_text("((A,B)0.5,C,D);((E,F)1.0,G,H);\n")
            reports = self.summaries(self.run_info(path))
        self.assertEqual([report["#"] for report in reports], ["1", "2"])
        self.assertAlmostEqual(float(reports[0]["MRP"]), 0.5)
        self.assertAlmostEqual(float(reports[1]["MRP"]), 0.99)
        self.assertAlmostEqual(float(reports[0]["MCS"]), -math.log(0.5))
        self.assertAlmostEqual(float(reports[1]["MCS"]), -math.log(0.01))

    def test_multiple_tree_standard_deviations_are_reported_per_tree(self):
        trees = (
            ("((A,B)0.2,(C,D)0.8,E);", [0.2, 0.8]),
            ("((A,B)20,(C,D)40,(E,F)60,G);", [0.2, 0.4, 0.6]),
            ("((A,B)1.0,C,D);", [0.99]),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            path.write_text("".join(tree for tree, _ in trees))
            reports = self.summaries(self.run_info(path))

        for label, transform in (
            ("sRP", lambda rp: rp),
            ("sCS", lambda rp: -math.log1p(-rp)),
        ):
            for report, (_, supports) in zip(reports, trees):
                values = [transform(rp) for rp in supports]
                mean = sum(values) / len(values)
                expected = math.sqrt(
                    sum((value - mean) ** 2 for value in values) / len(values)
                )
                with self.subTest(label=label, tree=report["#"]):
                    self.assertAlmostEqual(float(report[label]), expected)

    def test_tree_standard_deviations_are_na_only_for_unsupported_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trees.nwk"
            for trees in (
                "((A,B)0.2,(C,D)0.8,E);(A,B,C,D);",
                "(A,B,C,D);((A,B)0.2,(C,D)0.8,E);",
            ):
                with self.subTest(trees=trees):
                    path.write_text(trees)
                    reports = self.summaries(self.run_info(path))
                    supported = next(report for report in reports if report["support"] == "true")
                    unsupported = next(report for report in reports if report["support"] == "false")
                    self.assertNotEqual(supported["MRP"], "NA")
                    for column in ("MRP", "sRP", "MCS", "sCS"):
                        self.assertEqual(unsupported[column], "NA")


if __name__ == "__main__":
    unittest.main()
