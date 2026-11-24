AUTOMAKE_OPTIONS = subdir-objects
bin_PROGRAMS = get_newick_info
get_newick_info_SOURCES = src/get_newick_info.c
get_newick_info_CFLAGS = -std=c99 -O2 -Wall -Wextra -pedantic
get_newick_info_LDADD = -lm
EXTRA_DIST = CMakeLists.txt cmake/uninstall.cmake.in autogen.sh LICENSE \
	data/brlens.tree data/brlens_bp.tree data/brlens_bp.trees \
	data/brlens_pp.tree data/ultrametric.tree
