# get\_newick\_info

A command-line program for reporting statistics about tree files in [Newick
format](https://en.wikipedia.org/wiki/Newick_format).

## Download and Build

```sh
git clone https://github.com/nylander/get_newick_info.git
cd get_newick_info/src
make
```

## Usage examples

```sh
./get_newick_info ../data/brlens_pp.tree
./get_newick_info --taxa ../data/brlens_pp.tree
./get_newick_info --support ../data/brlens_pp.tree
./get_newick_info --tl ../data/brlens_pp.tree
./get_newick_info --tol 1e-7 ../data/ultrametric.tree
./get_newick_info --tol 1e-11 ../data/ultrametric.tree
```

Without an output option, the full summary is printed. `--taxa` (or `-t`)
prints sorted, unique taxon labels, `--support` prints support values, and
`--tl` prints the total branch length of each tree. These options print only
values, one per line, without headers or added indentation. These three output
options are mutually exclusive; use at most one per invocation.

Multiple input files are accepted:

```sh
./get_newick_info ../data/*.tree
```

With `--tl`, one value is printed per tree, in input order for each file.
The default summary instead reports the sum of tree lengths for each file. **TODO:** change to average if multi-tree file?

## Reported fields

* Number of trees
* Number of unique tip labels across all trees
* Whether branch lengths are present
* Rooted/unrooted status
* Whether all trees are dichotomous
* Whether comments are present
* Whether all trees are ultrametric
* Tree height when all trees are ultrametric and agree on height
* Total tree length (sum of all branch lengths)
* Mean resampling proportion (MRP)
* Mean character support (MCS)
* Optional sorted taxon labels
* Optional numeric internal-node support values
* Sum of support values

## Important conventions

### Rootedness

**TODO**: Check and try examples of "rooted" trees. Compare with Ape in R.

Plain Newick does not universally encode rootedness. This program therefore uses:

1. `[&R]` => rooted
2. `[&U]` => unrooted
3. Otherwise, a root with exactly two children is considered rooted.

This is a heuristic for files without explicit root metadata.

### Support

A numeric label immediately following an internal `)` is treated as a support
value. For example:

```text
((A:1,B:1)95:2,C:3)100;
```

recognizes `95` and `100` as support values.

Non-numeric internal labels are retained as node labels but are not included in
the support statistics.

MRP and MCS are calculated from all numeric support labels found in a tree,
using `n - 3` as the denominator regardless of how many support labels are
present, where `n` is the tip count. Support values whose maximum exceeds 1 are
treated as percentages and divided by 100. Values are limited to the range
0–1, and an RP of 1 is reduced to 0.99 before calculation. The summary reports
`NA` when no numeric support labels are found or the tree has three or fewer
tips. For files with multiple trees, MRP and MCS are the arithmetic means of
the per-tree metrics; if any tree has no usable support metrics, both
file-level metrics are `NA`.

### Ultrametricity

A tree is considered ultrametric when:

* every root-to-tip path has a branch length on every edge, and
* all root-to-tip distances agree within the specified relative/absolute
  tolerance.

The default tolerance is `1e-8`.

For a file containing multiple trees, the file-level `Is the tree ultrametric`
field is true only when every tree is ultrametric and their heights agree within
the tolerance. Otherwise `Tree height: NA` is reported.

### Comments

Any bracketed comment is detected. This includes ordinary Newick comments,
NHX annotations such as:

```text
[&&NHX:S=human:E=Homo_sapiens]
```

and BEAST/Enhanced-Newick-style comments such as:

```text
[&R]
[&U]
[&height=1.2]
```

## Examples

Input:

```text
((A:1,B:1)95:2,(C:1,D:1)90:2)100;
```

To print general information:

```sh
./get_newick_info example.nwk
```

Output:

```text
File: example.nwk
Number of trees: 1
Number of unique tip labels: 4
Branch lengths present: true
Rooted tree: true
Dichotomous: true
Comments present (NHX or Enhanced Newick format): false
Is the tree ultrametric: true
Tree height: 3
Tree length (TL): 8
Mean resampling proportion (MRP): NA
Mean character support (MCS): NA
Sum of support values: 285
```

To print taxon or tip labels:

```sh
./get_newick_info --taxa example.nwk
```

Output:

```text
A
B
C
D
```

To print support values:

```sh
./get_newick_info --support example.nwk
```

Output:

```text
95
90
100
```

To print only tree length (sum of all branch lengths):

```sh
./get_newick_info --tl example.nwk
```

```text
8
```

## Tests

Run the output regression tests (requires Python 3):

```sh
make test
```

## Limitations / design notes

For very deeply nested trees, the current implementation uses recursion while
constructing the tree. Typical phylogenetic trees are fine, but pathological
nesting can exhaust the C call stack.

## License

MIT License

## Download

<https://github.com/nylander/get_newick_info>
