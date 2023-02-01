import sys
import argparse
from gdscore.iomod import readintervalfile, savegsi, readgtreefile
from gdscore.rme import genPaszekGoreckiIntervals, genLCAIntervals, genGMSIntervals, genFHSIntervals, rme
from gdscore.rme_fhs import ppscores, merfellows
from gdscore.rec import rec
from gdscore.treeop import Tree, str2tree

MOD_LCA = 1                  # LCA model
MOD_GUIGO = 2                # GMS model Guigo et al.
MOD_PASZEKGORECKI = 3        # PG model
MOD_FELLOWS = 4              # FHS model Fellows et al.


def main():
    verbose = ""
    st = None
    gtrees = []
    gtreefile = None
    streefile = None
    outputfile = None
    # parse command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("-L", "--LCAmodel", action="store_const", const=1, default=0,
                        help="lca-mapping model of valid mappings (least flexible)")
    parser.add_argument("-G", "--GMSmodel", action="store_const", const=2, default=0,
                        help="GMS model of valid mappings")
    parser.add_argument("-P", "--PGmodel", action="store_const", const=3, default=0,
                        help="PG model of valid mappings (used by default if no flag is given)")
    parser.add_argument("-F", "--FHSmodel", action="store_const", const=4, default=0,
                        help="FHS (unrestricted) model of valid mappings")
    parser.add_argument("-i", "--setintervalsfromfile", action="store_const", const=1, default=0,
                        help="user-defined model of valid mappings; intervals are not generated from lca-mappings; "
                             "instead intervals are loaded from the input file")
    parser.add_argument("input_file", type=str, nargs="?", default=None, help="input file")
    parser.add_argument("-k", "--printtreewithscores", action="store_const", const=1, default=0,
                        help="print total computed score and a species tree with partial scores at corresponding nodes")
    parser.add_argument("-s", "--stree", type=str, nargs="?", help="-s TREE to define input species treee")
    parser.add_argument("-g", "--gtree", type=str, nargs="?", help="-g TREE to define input gene tree")
    parser.add_argument("-p", "--ptree", type=str, nargs="?",
                        help="-p 'TREE TREE' to give pair gene tree, species tree")
    parser.add_argument("-m", "--gtfile", type=str, nargs="?", help="-m FILE defines a set of gene trees in a file")
    parser.add_argument("-n", "--stfile", type=str, nargs="?", help="-n FILE defines a species tree from a file")
    parser.add_argument("-O", "--outfile", type=str, nargs="?",
                        help="-O FILE write detailed output to user defined output file")
    parser.add_argument("-j", "--selection", type=str, nargs="?",
                        help=" -j NUMTREE,LEN - skip first NUMTREE trees and take next LEN trees;"
                             " missing LEN=all trees")
    parser.add_argument("-d", "--delnodes", type=str, nargs="?",
                        help=" -d STRING      - list of node numbers to delete from SPEC and convert into duplications")
    parser.add_argument("-v", "--verbose", type=str, nargs="?",
                        help="(default) print ME score only; debug 0 - input/output; 1 - upper, lower bounds for "
                             "ME score; 2 - detailed Alg.2 debug; 3 - Alg.3 debug; 4 - skipped scores;"
                             "5 - print nodes from SPEC; 6 - more Alg.3 details; "
                             "7 - print size of SPEC and stops if > 22; 8 - print node number for species tree and "
                             "speciations from SPEC; where SPEC is a set of all lca-speciation nodes (from input gene "
                             "trees) which are ancestor of some duplication.")
    parser.add_argument("-t", "--showgtheights", action="store_const", const=1, default=0,
                        help="for each gene tree print its height")
    parser.add_argument("-e", "--computeECscore", action="store_const", const=1, default=0,
                        help="to compute EC score instead of ME score")

    args = parser.parse_args()
    printstreewithscores = args.printtreewithscores  # by default we do not print tree with scores
    model = args.LCAmodel + args.GMSmodel + args.PGmodel + args.FHSmodel  # all not choosen flags equal 0
    if model == 0:
        model = 3  # if no model is set, the default model is PG
    if args.setintervalsfromfile:
        model = 0  # model is set to custom model loaded from file
    if args.stree:
        st = Tree(str2tree(args.stree))  # species tree given from input
    if args.gtree:
        gtrees.append(Tree(str2tree(args.gtree)))  # gene tree given from input
    if args.ptree:  # pair of gene tree species tree given
        gtree, stree = args.ptree.strip().split(" ")
        gtrees.append(Tree(str2tree(gtree)))
        st = Tree(str2tree(stree))
    if args.gtfile:
        gtreefile = args.gtfile
    if args.stfile:
        streefile = args.stfile
    if args.outfile:
        outputfile = args.outfile
    firstgtrees = -1
    numtrees = -1
    if args.selection:
        if "," in args.selection:
            firstgtrees, numtrees = int(args.selection.split(",")[0]), int(args.selection.split(",")[1])
        else:
            firstgtrees, numtrees = int(args.selection), -1
    delnodes = ""
    if args.delnodes:
        delnodes = args.delnodes
    if args.verbose:
        verbose = args.verbose

    if gtreefile:  # gene trees are given in a separate file, -m option
        gtrees.extend(readgtreefile(gtreefile))
    if streefile:
        treelist = readgtreefile(streefile)
        if treelist:
            st = treelist[0]
    if args.input_file:
        gtrees, st = readintervalfile(args.input_file)  # gene trees and species tree obtained from input file
    if args.showgtheights:
        for gt in gtrees:
            print(gt.height(), gt)

    if not st or not gtrees:  # missing tree
        print("Both trees have to be defined")  # ex usage of only -s option // or -m without -s ect
        sys.exit(3)

    runmer = 1  # default run

    # for the selection option
    if firstgtrees >= 0:
        if numtrees > 0:
            gtrees = gtrees[firstgtrees:firstgtrees + numtrees]
        else:
            gtrees = gtrees[firstgtrees:]

    # sets the lca-mappings for gene tree nodes (skipped for custom set intervals)
    if model:
        for gt in gtrees:
            gt.set_lca_mapping(st)  # changed to embretnet variant

    # generates intervals according to chosen model
    for gt in gtrees:
        if model == MOD_PASZEKGORECKI:
            genPaszekGoreckiIntervals(gt, st)
        elif model == MOD_GUIGO:
            genGMSIntervals(gt, st)
        elif model == MOD_LCA:
            genLCAIntervals(gt, st)

    if args.computeECscore:
        if model == MOD_FELLOWS:  # creates intervals for all nodes
            for gt in gtrees:
                genFHSIntervals(gt, st)
        ec = rec(gtrees, st)
        print("EC score", ec)
        sys.exit(0)

    if model == MOD_FELLOWS:
        me = merfellows(gtrees, st, verbose, printstreewithscores)
        print("ME score", me)
        sys.exit(0)

    if outputfile:
        savegsi(gtrees, st, outputfile)

    if verbose.find("3") != -1:
        for i, gt in enumerate(gtrees):
            for g in gt.root.nodes():
                if g.interval:
                    print("Interval for %d in tree%d:" % (g.num, i), g, g.interval)

    if runmer:
        me = rme(gtrees, st, verbose)
        if printstreewithscores:
            print("&s", ppscores(st))
        else:
            print("ME score", me)

    return 0


if __name__ == "__main__":
    main()
