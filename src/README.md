# newick-info

A small standalone C99 command-line program for reporting statistics about
Newick tree files.

## Build

```sh
make
```

or:

```sh
cc -std=c99 -O2 -Wall -Wextra -pedantic -o newick-info newick-info.c -lm
```

## Usage

```sh
./newick-info trees.nwk
./newick-info --taxa trees.nwk
./newick-info --support trees.nwk
./newick-info --taxa --support trees.nwk
./newick-info --tol 1e-7 trees.nwk
```

Multiple input files are accepted:

```sh
./newick-info a.nwk b.nwk c.nwk
```

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
./newick-info --taxa --support example.nwk
```

Expected report is approximately:

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
Tree length, TL: 10
Taxon labels:
  A
  B
  C
  D
Support values:
  95
  90
  100
Sum of support values: 285
```

## Limitations / design notes

This is intentionally a small parser rather than a full NEXUS parser. It reads
Newick trees from the supplied files; it does not parse NEXUS `BEGIN TREES`
blocks or translate tables.

For very deeply nested trees, the current implementation uses recursion while
constructing the tree. Typical phylogenetic trees are fine, but pathological
nesting can exhaust the C call stack.
