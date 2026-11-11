# get_newick_info

A small standalone C99 command-line program for reporting statistics about
Newick tree files.

## Build

```sh
make
```

or:

```sh
cc -std=c99 -O2 -Wall -Wextra -pedantic -o get_newick_info get_newick_info.c -lm
```

## Usage

```sh
./get_newick_info trees.nwk
./get_newick_info --taxa trees.nwk
./get_newick_info --support trees.nwk
./get_newick_info --tl trees.nwk
./get_newick_info --taxa --support trees.nwk
./get_newick_info --tol 1e-7 trees.nwk
```

Without an output option, the full summary is printed. `--taxa` (or `-t`)
prints sorted, unique taxon labels, `--support` prints support values, and
`--tl` prints the total branch length of each tree. These options print only
values, one per line, without headers or added indentation. When combined,
values appear in that order: taxa, support, then tree length.

Multiple input files are accepted:

```sh
./get_newick_info a.nwk b.nwk c.nwk
```

With `--tl`, one value is printed per tree, in input order for each file.
The default summary instead reports the sum of tree lengths for each file.

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
* Optional sorted taxon labels
* Optional numeric internal-node support values
* Sum of support values

## Important conventions

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

## Example

Input:

```text
((A:1,B:1)95:2,(C:1,D:1)90:2)100;
```

Command:

```sh
./get_newick_info --taxa --support example.nwk
```

Output:

```text
A
B
C
D
95
90
100
```

To print only tree length:

```sh
./get_newick_info --tl example.nwk
```

```text
10
```

## Tests

Run the output regression tests (requires Python 3):

```sh
make test
```

## Limitations / design notes

This is intentionally a small parser rather than a full NEXUS parser. It reads
Newick trees from the supplied files; it does not parse NEXUS `BEGIN TREES`
blocks or translate tables.

For very deeply nested trees, the current implementation uses recursion while
constructing the tree. Typical phylogenetic trees are fine, but pathological
nesting can exhaust the C call stack.
