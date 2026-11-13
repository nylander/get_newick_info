# get\_newick\_info

A command-line program for reporting statistics about tree files in [Newick
format](https://en.wikipedia.org/wiki/Newick_format).

![](img/img.jpg)

---

## Reported fields

| Column | Label       | Property |
| ------ | ----------- | -------- |
| 1.     | **#**       |  Position of this tree in the input file |
| 2.     | **Ntips**   |  Number of tips in this tree |
| 3.     | **blens**   |  Are branch lengths present? |
| 4.     | **support** |  Do branches have support values? |
| 5.     | **root**    |  Is the tree rooted? |
| 6.     | **dicho**   |  Is the tree dichotomous? |
| 7.     | **comm**    |  Are there comments present? |
| 8.     | **ultra**   |  Is the tree ultrametric? |
| 9.     | **height**  |  Tree height (when ultrametric) |
| 10.    | **TL**      |  Total length of this tree's branches |
| 11.    | **sumSup**  |  Sum of support values |
| 12.    | **MRP**     |  Mean resampling proportion |
| 13.    | **sRP**     |  Standard deviation of resampling proportion |
| 14.    | **MCS**     |  Mean character support |
| 15.    | **sCS**     |  Standard deviation of character support |
| 16.    | **Ntrees**  |  Total number of trees in the input file |
| 17.    | **file**    |  Input file name |

## Usage and Options

    Usage: ./get_newick_info [options] file1.nwk [file2.nwk ...]

    Options:
      -t,--taxa        print sorted taxon labels
      -s, --support    print internal-node support values
      -b, --brlens     print branch lengths
      --tl             print tree length
      --tol VALUE      ultrametric comparison tolerance (default 1e-8)
      --lp VALUE       correction value for calculating MRP and MCS (default 0.01)
      -n, --noverbose  do not print header in output
      -p, --fullpath   print realpath to input file
      -v, --version    print version
      -h, --help       show this help

    Use at most one of --taxa, --support, --brlens and --tl.

## Download and Build

Note: requires a C compiler (and python3 for testing)

```sh
git clone https://github.com/nylander/get_newick_info.git
cd get_newick_info/src
make
```

## Usage examples

```sh
./get_newick_info ../data/brlens_pp.tree
./get_newick_info -n ../data/brlens_pp.tree
./get_newick_info --taxa ../data/brlens_pp.tree
./get_newick_info --support ../data/brlens_pp.tree
./get_newick_info --brlens ../data/brlens_pp.tree
./get_newick_info --tl ../data/brlens_pp.tree
./get_newick_info --tol 1e-7 ../data/ultrametric.tree
./get_newick_info --tol 1e-11 ../data/ultrametric.tree
```

Without an output option, the full summary is printed. `--taxa` (or `-t`)
prints sorted, unique taxon labels, `--support` prints support values,
`--brlens` prints the branch lengths, and `--tl` prints the total branch length
of each tree. These options print only values, one per line, without headers or
added indentation. These three output options are mutually exclusive; use at
most one per invocation.

Multiple input files are accepted. Example:

```sh
./get_newick_info -n ../data/*.tree | column -t | less -S
```

Multi-tree files are accepted. Each tree is reported on its own row, with its
position in the file in the first column and the file's total tree count in the
second-to-last column. Example:

```sh
./get_newick_info -n ../data/*.trees | column -t | less -S
```


## Important conventions

Note: some features are still untested, see [TODO](#todo) below.

### Rootedness

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

### MCS and MRP

MRP and MCS are calculated from all numeric support labels found in a tree,
using `n - 3` as the denominator regardless of how many support labels are
present, where `n` is the tip count. Support values whose maximum exceeds 1 are
treated as percentages and divided by 100. Values are limited to the range
0--1, and an RP of 1 is reduced to $1-lp$ before calculation (to avoid division by 0).
The value of `lp` can be set by the `--lp` option and default is 0.01.

Each tree is reported separately. MRP and MCS are `NA` when no numeric support
labels are found or the tree has three or fewer tips.

The per-tree standard deviations `s_RP` and `s_CS` are population standard
deviations over the available support labels, dividing by the number of labels
rather than `n - 3`. They use the same normalized RP values as MRP and the
transformed values `CS = -log(1 - RP)` as MCS. Both standard deviations are
`NA` when the individual tree has no usable support metrics.

### Ultrametricity

A tree is considered ultrametric when:

* every root-to-tip path has a branch length on every edge, and
* all root-to-tip distances agree within the specified relative/absolute
  tolerance.

The default tolerance is `1e-8`.

Ultrametricity and height are reported separately for each tree.

### Comments

(See also [TODO](#todo) below)

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

MIT [License](LICENSE)

## Download

<https://github.com/nylander/get_newick_info>

## TODO

- [ ] Test with BEAST-extended files
- [ ] Test with almost-ultrametric files
- [ ] Test with rooted/unrooted trees
- [ ] Test with non-dichotomous trees
- [ ] Manually check MRP and MCS calculations on multi-tree files
- [ ] Add reference to MRP and MCS

