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
            self.assertIn("Tree length, TL: 10\n", self.run_info(path))
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
            f"Tree length, TL: {TL}"
            "Sum of support values: 274\n\n",
        )


if __name__ == "__main__":
    unittest.main()
