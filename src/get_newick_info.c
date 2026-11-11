
/*
 * newick-info.c
 *
 * C99 program for reporting statistics about Newick trees.
 *
 * Build:
 *   cc -std=c99 -O2 -Wall -Wextra -pedantic -o newick-info newick-info.c -lm
 *
 * Usage:
 *   newick-info [options] file1.nwk [file2.nwk ...]
 *
 * Options:
 *   --taxa, -t       print sorted taxon labels
 *   --support       print recognized numeric internal-node support values
 *   --tl            print tree length for each tree
 *   --tol VALUE      tolerance for ultrametric comparisons (default 1e-8)
 *   --help, -h       show help
 *
 * Support convention:
 *   A numeric label on an internal node is treated as a support value.
 *   Non-numeric internal labels are not treated as support.
 *
 * Rootedness convention:
 *   Newick does not, by itself, universally encode rootedness. This program
 *   reports rooted=true when an explicit BEAST-style [&R] root annotation is
 *   present, rooted=false when [&U] is present, otherwise it uses the common
 *   structural heuristic that a root with exactly two children is rooted.
 *
 * Comments:
 *   Any bracketed Newick comment is counted. This includes NHX and
 *   BEAST/Enhanced-Newick annotations.
 */

#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Node Node;

struct Node {
    char *label;
    double length;
    bool has_length;
    bool is_leaf;
    bool is_support;
    double support;
    bool has_comment;
    int child_count;
    Node **children;
};

typedef struct {
    const char *s;
    size_t pos;
    size_t len;
    int error;
    char error_msg[256];
} Parser;

typedef struct {
    size_t trees;
    size_t tips;
    size_t internal_nodes;
    size_t branch_count;
    size_t branch_lengths;
    double tree_length;
    double height;
    bool ultrametric;
    bool rooted;
    bool explicit_rootedness;
    bool explicit_unrootedness;
    bool comments;
    bool dichotomous;
    bool have_support;
    double support_sum;
    char **taxa;
    size_t ntaxa;
    size_t taxa_cap;
    double *supports;
    size_t nsupports;
    size_t support_cap;
} FileStats;

static void die_oom(void) {
    fprintf(stderr, "fatal: out of memory\n");
    exit(EXIT_FAILURE);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p && n) die_oom();
    return p;
}

static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s);
    if (!p && n && s) die_oom();
    return p;
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n);
    if (!q && n) die_oom();
    return q;
}

static char *xstrdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *p = xmalloc(n);
    memcpy(p, s, n);
    return p;
}

static void parser_error(Parser *p, const char *msg) {
    if (!p->error) {
        p->error = 1;
        snprintf(p->error_msg, sizeof p->error_msg,
                 "at character %zu: %s", p->pos + 1, msg);
    }
}

static void skip_ws(Parser *p) {
    while (p->pos < p->len && isspace((unsigned char)p->s[p->pos]))
        p->pos++;
}

static bool peek(Parser *p, char c) {
    skip_ws(p);
    return p->pos < p->len && p->s[p->pos] == c;
}

static bool consume(Parser *p, char c) {
    skip_ws(p);
    if (p->pos < p->len && p->s[p->pos] == c) {
        p->pos++;
        return true;
    }
    return false;
}

static char *read_quoted(Parser *p) {
    char quote = p->s[p->pos++];
    size_t cap = 64, n = 0;
    char *out = xmalloc(cap);

    while (p->pos < p->len) {
        char c = p->s[p->pos++];

        /* Newick/NEXUS convention: doubled quote escapes a quote. */
        if (c == quote && p->pos < p->len && p->s[p->pos] == quote) {
            p->pos++;
        } else if (c == quote) {
            out[n] = '\0';
            return out;
        }

        if (n + 2 > cap) {
            cap *= 2;
            out = xrealloc(out, cap);
        }
        out[n++] = c;
    }

    free(out);
    parser_error(p, "unterminated quoted label");
    return NULL;
}

static bool is_label_delim(char c) {
    return isspace((unsigned char)c) || c == '(' || c == ')' ||
           c == ',' || c == ':' || c == ';' || c == '[';
}

static char *read_label(Parser *p) {
    skip_ws(p);
    if (p->pos >= p->len) return NULL;

    if (p->s[p->pos] == '\'' || p->s[p->pos] == '"')
        return read_quoted(p);

    size_t start = p->pos;
    while (p->pos < p->len && !is_label_delim(p->s[p->pos]))
        p->pos++;

    if (p->pos == start) return NULL;

    size_t n = p->pos - start;
    char *label = xmalloc(n + 1);
    memcpy(label, p->s + start, n);
    label[n] = '\0';
    return label;
}

static bool parse_number(Parser *p, double *value) {
    skip_ws(p);
    const char *start = p->s + p->pos;
    char *end = NULL;

    errno = 0;
    double x = strtod(start, &end);

    if (end == start) return false;
    if (errno == ERANGE || !isfinite(x)) {
        parser_error(p, "invalid or out-of-range numeric value");
        return false;
    }

    p->pos += (size_t)(end - start);
    *value = x;
    return true;
}

static bool parse_support_number(const char *s, double *v) {
    if (!s || !*s) return false;

    char *end = NULL;
    errno = 0;
    double x = strtod(s, &end);
    if (end == s || *end != '\0' || errno == ERANGE || !isfinite(x))
        return false;

    *v = x;
    return true;
}

static void node_add_child(Node *parent, Node *child) {
    parent->children = xrealloc(parent->children,
                                (size_t)(parent->child_count + 1) * sizeof *parent->children);
    parent->children[parent->child_count++] = child;
}

static Node *node_new(bool leaf) {
    Node *n = xcalloc(1, sizeof *n);
    n->length = 0.0;
    n->is_leaf = leaf;
    return n;
}

static void node_free(Node *n) {
    if (!n) return;
    for (int i = 0; i < n->child_count; i++)
        node_free(n->children[i]);
    free(n->children);
    free(n->label);
    free(n);
}

static void consume_comment(Parser *p, FileStats *fs) {
    skip_ws(p);
    if (!peek(p, '[')) return;

    p->pos++; /* '[' */
    fs->comments = true;

    size_t start = p->pos;
    int depth = 1;

    while (p->pos < p->len && depth > 0) {
        char c = p->s[p->pos++];
        if (c == '[') depth++;
        else if (c == ']') depth--;
    }

    if (depth != 0) {
        parser_error(p, "unterminated bracket comment");
        return;
    }

    size_t n = (p->pos - 1) - start;
    char *comment = xmalloc(n + 1);
    memcpy(comment, p->s + start, n);
    comment[n] = '\0';

    /* Detect explicit BEAST-style root flags in comments. */
    for (size_t i = 0; i + 2 < n; i++) {
        if (comment[i] == '&' && comment[i + 1] == 'R') {
            fs->rooted = true;
            fs->explicit_rootedness = true;
        }
        if (comment[i] == '&' && comment[i + 1] == 'U') {
            fs->rooted = false;
            fs->explicit_unrootedness = true;
        }
    }

    free(comment);
}

static void consume_comments(Parser *p, FileStats *fs) {
    for (;;) {
        skip_ws(p);
        if (p->pos < p->len && p->s[p->pos] == '[')
            consume_comment(p, fs);
        else
            break;
    }
}

static Node *parse_subtree(Parser *p, FileStats *fs);

static Node *parse_internal(Parser *p, FileStats *fs) {
    Node *n = node_new(false);

    if (!consume(p, '(')) {
        node_free(n);
        parser_error(p, "expected '('");
        return NULL;
    }

    for (;;) {
        Node *child = parse_subtree(p, fs);
        if (!child) {
            node_free(n);
            return NULL;
        }
        node_add_child(n, child);

        if (consume(p, ','))
            continue;

        if (consume(p, ')'))
            break;

        node_free(n);
        parser_error(p, "expected ',' or ')'");
        return NULL;
    }

    consume_comments(p, fs);

    char *label = read_label(p);
    if (label) {
        n->label = label;
        double support;
        if (parse_support_number(label, &support)) {
            n->is_support = true;
            n->support = support;
            fs->have_support = true;
            fs->support_sum += support;

            if (fs->nsupports == fs->support_cap) {
                fs->support_cap = fs->support_cap ? fs->support_cap * 2 : 16;
                fs->supports = xrealloc(fs->supports,
                                         fs->support_cap * sizeof *fs->supports);
            }
            fs->supports[fs->nsupports++] = support;
        }
    }

    consume_comments(p, fs);

    if (consume(p, ':')) {
        double length;
        if (!parse_number(p, &length)) {
            node_free(n);
            parser_error(p, "expected branch length after ':'");
            return NULL;
        }
        n->length = length;
        n->has_length = true;
        fs->branch_lengths++;
        fs->branch_count++;
        fs->tree_length += length;
    }

    consume_comments(p, fs);
    return n;
}

static Node *parse_leaf(Parser *p, FileStats *fs) {
    Node *n = node_new(true);

    char *label = read_label(p);
    if (!label) {
        node_free(n);
        parser_error(p, "expected a taxon label");
        return NULL;
    }

    n->label = label;
    fs->tips++;

    if (fs->ntaxa == fs->taxa_cap) {
        fs->taxa_cap = fs->taxa_cap ? fs->taxa_cap * 2 : 32;
        fs->taxa = xrealloc(fs->taxa, fs->taxa_cap * sizeof *fs->taxa);
    }
    fs->taxa[fs->ntaxa++] = xstrdup(label);

    consume_comments(p, fs);

    if (consume(p, ':')) {
        double length;
        if (!parse_number(p, &length)) {
            node_free(n);
            parser_error(p, "expected branch length after ':'");
            return NULL;
        }
        n->length = length;
        n->has_length = true;
        fs->branch_lengths++;
        fs->branch_count++;
        fs->tree_length += length;
    }

    consume_comments(p, fs);
    return n;
}

static Node *parse_subtree(Parser *p, FileStats *fs) {
    skip_ws(p);
    if (p->pos >= p->len) {
        parser_error(p, "unexpected end of input");
        return NULL;
    }

    if (p->s[p->pos] == '(')
        return parse_internal(p, fs);

    return parse_leaf(p, fs);
}

static void reset_per_tree_stats(FileStats *fs) {
    fs->tips = 0;
    fs->internal_nodes = 0;
    fs->branch_count = 0;
    fs->branch_lengths = 0;
    fs->tree_length = 0.0;
    fs->height = 0.0;
    fs->ultrametric = true;
    fs->dichotomous = true;
    fs->have_support = false;
    fs->support_sum = 0.0;

    free(fs->taxa);
    fs->taxa = NULL;
    fs->ntaxa = fs->taxa_cap = 0;

    free(fs->supports);
    fs->supports = NULL;
    fs->nsupports = fs->support_cap = 0;
}

static bool all_branch_lengths_present(const Node *n) {
    if (n->is_leaf) return n->has_length;
    for (int i = 0; i < n->child_count; i++) {
        if (!n->children[i]->has_length) return false;
        if (!all_branch_lengths_present(n->children[i])) return false;
    }
    return true;
}

static double root_to_tip(const Node *n, double accumulated,
                          double *first_height, bool *first_set,
                          bool *ultrametric, double tol) {
    double here = accumulated;
    if (!n->is_leaf && n->has_length)
        here += n->length;

    if (n->is_leaf) {
        if (!n->has_length) {
            *ultrametric = false;
            return NAN;
        }
        here += n->length;

        if (!*first_set) {
            *first_height = here;
            *first_set = true;
        } else {
            double scale = fmax(1.0, fmax(fabs(*first_height), fabs(here)));
            if (fabs(*first_height - here) > tol * scale)
                *ultrametric = false;
        }
        return here;
    }

    double result = NAN;
    for (int i = 0; i < n->child_count; i++) {
        double x = root_to_tip(n->children[i], here, first_height,
                               first_set, ultrametric, tol);
        if (!isnan(x)) result = x;
    }
    return result;
}

/*
 * The root itself normally has no branch length in a rooted Newick tree.
 * Therefore root_to_tip starts at zero and adds branch lengths on children.
 */
static void analyze_tree(Node *root, FileStats *fs, double tol) {
    fs->dichotomous = true;

    /* Root structural heuristic, unless explicit [&R]/[&U] was found. */
    if (!fs->explicit_rootedness && !fs->explicit_unrootedness)
        fs->rooted = (root->child_count == 2);

    /* Internal-node count and dichotomy. */
    Node **stack = NULL;
    size_t nstack = 0, cap = 0;

    cap = 64;
    stack = xmalloc(cap * sizeof *stack);
    stack[nstack++] = root;

    while (nstack) {
        Node *n = stack[--nstack];

        if (!n->is_leaf) {
            fs->internal_nodes++;
            if (n->child_count != 2)
                fs->dichotomous = false;

            for (int i = 0; i < n->child_count; i++) {
                if (nstack == cap) {
                    cap *= 2;
                    stack = xrealloc(stack, cap * sizeof *stack);
                }
                stack[nstack++] = n->children[i];
            }
        }
    }
    free(stack);

    bool all_lengths = all_branch_lengths_present(root);
    (void)all_lengths;

    if (fs->branch_lengths == 0) {
        fs->ultrametric = false;
        fs->height = NAN;
    } else {
        bool first_set = false;
        double first_height = NAN;
        bool ultra = true;
        (void)root_to_tip(root, 0.0, &first_height, &first_set, &ultra, tol);

        fs->ultrametric = ultra && first_set && all_branch_lengths_present(root);
        fs->height = fs->ultrametric ? first_height : NAN;
    }
}

static int cmp_str(const void *a, const void *b) {
    const char * const *sa = a;
    const char * const *sb = b;
    return strcmp(*sa, *sb);
}

static size_t count_unique_labels(char **labels, size_t n) {
    if (n == 0) return 0;
    qsort(labels, n, sizeof *labels, cmp_str);
    size_t unique = 1;
    for (size_t i = 1; i < n; i++) {
        if (strcmp(labels[i - 1], labels[i]) != 0)
            unique++;
    }
    return unique;
}

static void print_bool(bool b) {
    puts(b ? "true" : "false");
}

static void print_file_report(const char *filename, const FileStats *fs,
                              bool print_taxa, bool print_support, bool print_tl,
                              const double *tree_lengths) {
    if (!print_taxa && !print_support && !print_tl) {
        printf("File: %s\n", filename);
        printf("Number of trees: %zu\n", fs->trees);
        printf("Number of unique tip labels: %zu\n", fs->tips);
        printf("Branch lengths present: ");
        print_bool(fs->branch_lengths > 0);
        printf("Rooted tree: ");
        print_bool(fs->rooted);
        printf("Dichotomous: ");
        print_bool(fs->dichotomous);
        printf("Comments present (NHX or Enhanced Newick format): ");
        print_bool(fs->comments);

        /*
         * For a multi-tree file, the summary is meaningful only if every tree
         * is ultrametric. The caller aggregates the per-tree result below.
         */
        printf("Is the tree ultrametric: ");
        print_bool(fs->ultrametric);

        if (fs->ultrametric)
            printf("Tree height: %.17g\n", fs->height);
        else
            printf("Tree height: NA\n");

        printf("Tree length, TL: %.17g\n", fs->tree_length);
    }

    if (print_taxa) {
        char **copy = xmalloc(fs->ntaxa * sizeof *copy);
        for (size_t i = 0; i < fs->ntaxa; i++) copy[i] = fs->taxa[i];
        qsort(copy, fs->ntaxa, sizeof *copy, cmp_str);

        for (size_t i = 0; i < fs->ntaxa; i++) {
            if (i == 0 || strcmp(copy[i], copy[i - 1]) != 0)
                printf("%s\n", copy[i]);
        }

        free(copy);
    }

    if (print_support) {
        for (size_t i = 0; i < fs->nsupports; i++)
            printf("%.17g\n", fs->supports[i]);
    }
    if (print_tl) {
        for (size_t i = 0; i < fs->trees; i++)
            printf("%.17g\n", tree_lengths[i]);
    }

    if (!print_taxa && !print_support && !print_tl) {
        printf("Sum of support values: %.17g\n", fs->support_sum);
        putchar('\n');
    }
}

static char *read_entire_file(const char *filename, size_t *len_out) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        return NULL;
    }

    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        fprintf(stderr, "%s: cannot seek\n", filename);
        return NULL;
    }

    long n = ftell(fp);
    if (n < 0) {
        fclose(fp);
        fprintf(stderr, "%s: cannot determine file size\n", filename);
        return NULL;
    }
    rewind(fp);

    char *buf = xmalloc((size_t)n + 1);
    size_t got = fread(buf, 1, (size_t)n, fp);
    if (ferror(fp)) {
        fprintf(stderr, "%s: read error\n", filename);
        free(buf);
        fclose(fp);
        return NULL;
    }
    fclose(fp);

    buf[got] = '\0';
    *len_out = got;
    return buf;
}

static void free_stats(FileStats *fs) {
    for (size_t i = 0; i < fs->ntaxa; i++)
        free(fs->taxa[i]);
    free(fs->taxa);
    free(fs->supports);
}

static bool parse_file(const char *filename, bool print_taxa, bool print_support,
                       bool print_tl, double tol) {
    size_t len;
    char *text = read_entire_file(filename, &len);
    if (!text) return false;

    Parser p = { text, 0, len, 0, {0} };
    FileStats fs = {0};

    /*
     * We accumulate file-wide totals, while retaining the most recent
     * tree's detailed values. For multiple trees, we separately track
     * whether all trees are ultrametric and whether their heights agree.
     */
    size_t total_branch_lengths = 0;
    double total_tree_length = 0.0;
    double total_support_sum = 0.0;
    size_t total_supports = 0;
    bool any_comments = false;
    bool all_dichotomous = true;
    bool all_ultrametric = true;
    bool first_height_set = false;
    double common_height = NAN;
    bool heights_agree = true;
    bool rooted_value_set = false;
    bool rooted_value = false;

    char **all_taxa = NULL;
    size_t all_ntaxa = 0, all_taxa_cap = 0;
    double *all_supports = NULL;
    size_t all_nsupports = 0, all_support_cap = 0;
    double *tree_lengths = NULL;
    size_t tree_lengths_cap = 0;

    while (1) {
        skip_ws(&p);
        if (p.pos >= p.len) break;

        reset_per_tree_stats(&fs);

        Node *root = parse_subtree(&p, &fs);
        if (!root) {
            fprintf(stderr, "%s: parse error %s\n", filename, p.error_msg);
            node_free(root);
            free(text);
            free_stats(&fs);
            free(all_taxa);
            free(all_supports);
            free(tree_lengths);
            return false;
        }

        consume_comments(&p, &fs);

        if (!consume(&p, ';')) {
            node_free(root);
            fprintf(stderr, "%s: parse error at character %zu: expected ';'\n",
                    filename, p.pos + 1);
            free(text);
            free_stats(&fs);
            free(all_taxa);
            free(all_supports);
            free(tree_lengths);
            return false;
        }

        analyze_tree(root, &fs, tol);

        if (print_tl) {
            if (fs.trees == tree_lengths_cap) {
                tree_lengths_cap = tree_lengths_cap ? tree_lengths_cap * 2 : 16;
                tree_lengths = xrealloc(tree_lengths,
                                        tree_lengths_cap * sizeof *tree_lengths);
            }
            tree_lengths[fs.trees] = fs.tree_length;
        }
        fs.trees++;
        total_branch_lengths += fs.branch_lengths;
        total_tree_length += fs.tree_length;
        total_support_sum += fs.support_sum;
        total_supports += fs.nsupports;
        any_comments = any_comments || fs.comments;
        all_dichotomous = all_dichotomous && fs.dichotomous;
        all_ultrametric = all_ultrametric && fs.ultrametric;

        if (fs.trees == 1) {
            rooted_value = fs.rooted;
            rooted_value_set = true;
        } else if (rooted_value_set && rooted_value != fs.rooted) {
            /*
             * A file containing mixed explicit/structural rootedness is
             * reported as false rather than inventing a third state.
             */
            rooted_value = false;
        }

        if (fs.ultrametric) {
            if (!first_height_set) {
                common_height = fs.height;
                first_height_set = true;
            } else {
                double scale = fmax(1.0, fmax(fabs(common_height), fabs(fs.height)));
                if (fabs(common_height - fs.height) > tol * scale)
                    heights_agree = false;
            }
        }

        for (size_t i = 0; i < fs.ntaxa; i++) {
            if (all_ntaxa == all_taxa_cap) {
                all_taxa_cap = all_taxa_cap ? all_taxa_cap * 2 : 64;
                all_taxa = xrealloc(all_taxa,
                                    all_taxa_cap * sizeof *all_taxa);
            }
            all_taxa[all_ntaxa++] = xstrdup(fs.taxa[i]);
        }

        if (print_support) {
            for (size_t i = 0; i < fs.nsupports; i++) {
                if (all_nsupports == all_support_cap) {
                    all_support_cap = all_support_cap ? all_support_cap * 2 : 64;
                    all_supports = xrealloc(all_supports,
                                             all_support_cap * sizeof *all_supports);
                }
                all_supports[all_nsupports++] = fs.supports[i];
            }
        }

        node_free(root);
    }

    /* Construct a file-level report from accumulated values. */
    FileStats report = {0};
    report.trees = fs.trees;
    report.tips = count_unique_labels(all_taxa, all_ntaxa);
    report.branch_lengths = total_branch_lengths;
    report.tree_length = total_tree_length;
    report.comments = any_comments;
    report.dichotomous = all_dichotomous;
    report.ultrametric = all_ultrametric && first_height_set && heights_agree;
    report.height = report.ultrametric ? common_height : NAN;
    report.rooted = rooted_value_set ? rooted_value : false;
    report.support_sum = total_support_sum;
    report.ntaxa = all_ntaxa;
    report.taxa = all_taxa;
    report.nsupports = all_nsupports;
    report.supports = all_supports;

    print_file_report(filename, &report, print_taxa, print_support, print_tl,
                      tree_lengths);

    free(text);
    free(tree_lengths);
    free_stats(&fs);
    free_stats(&report);
    return true;
}

static void usage(const char *prog) {
    printf("Usage: %s [options] file1.nwk [file2.nwk ...]\n", prog);
    puts("");
    puts("Options:");
    puts("  --taxa, -t       print sorted taxon labels");
    puts("  --support        print numeric internal-node support values");
    puts("  --tl             print tree length for each tree");
    puts("  --tol VALUE      ultrametric comparison tolerance (default 1e-8)");
    puts("  --help, -h       show this help");
    puts("");
    puts("Use at most one of --taxa, --support, and --tl.");
    puts("Numeric internal labels are interpreted as support values.");
}

int main(int argc, char **argv) {
    bool print_taxa = false;
    bool print_support = false;
    bool print_tl = false;
    double tol = 1e-8;

    int first_file = 1;

    while (first_file < argc) {
        const char *arg = argv[first_file];

        if (strcmp(arg, "--taxa") == 0 || strcmp(arg, "-t") == 0) {
            print_taxa = true;
            first_file++;
        } else if (strcmp(arg, "--support") == 0) {
            print_support = true;
            first_file++;
        } else if (strcmp(arg, "--tl") == 0) {
            print_tl = true;
            first_file++;
        } else if (strcmp(arg, "--tol") == 0) {
            if (first_file + 1 >= argc) {
                fprintf(stderr, "--tol requires a value\n");
                return EXIT_FAILURE;
            }
            char *end = NULL;
            errno = 0;
            tol = strtod(argv[first_file + 1], &end);
            if (errno || end == argv[first_file + 1] || *end != '\0' ||
                tol < 0.0 || !isfinite(tol)) {
                fprintf(stderr, "invalid --tol value: %s\n", argv[first_file + 1]);
                return EXIT_FAILURE;
            }
            first_file += 2;
        } else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0) {
            usage(argv[0]);
            return EXIT_SUCCESS;
        } else if (arg[0] == '-') {
            fprintf(stderr, "unknown option: %s\n", arg);
            usage(argv[0]);
            return EXIT_FAILURE;
        } else {
            break;
        }
    }

    if ((print_taxa ? 1 : 0) + (print_support ? 1 : 0) + (print_tl ? 1 : 0) > 1) {
        fprintf(stderr, "--taxa, --support, and --tl are mutually exclusive\n");
        return EXIT_FAILURE;
    }

    if (first_file >= argc) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    int failures = 0;
    for (int i = first_file; i < argc; i++) {
        if (!parse_file(argv[i], print_taxa, print_support, print_tl, tol))
            failures++;
    }

    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
