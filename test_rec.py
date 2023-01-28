import pytest
from iomod import readgtreefile, readintervalfile
from rme import genFHSIntervals, genLCAIntervals, genGMSIntervals, genPaszekGoreckiIntervals
from rec import rec

MOD_LCA = 1  # LCA model
MOD_GUIGO = 2  # GMS model Guigo et al.
MOD_PASZEKGORECKI = 3  # PG model
MOD_FELLOWS = 4  # FHS model Fellows et al.

# Test corresponds to published results in
# J. Paszek and P. Gorecki, “Genomic duplication problems for unrooted gene trees,”
# BMC Genomics, vol. 17, no. 1, pp. 165–175, 2016.
# The results are based on table 1 and extended
# Note here we report 4 instead of 5 for GMS
infiles = ["data/RME/in.txt", "data/RME/in1.txt", "data/RME/in2.txt"]
inrecfiles = ["data/REC/guigo/g.txt", "data/REC/guigo/s1.txt",
              "data/REC/treefam/g.txt", "data/REC/treefam/s.txt"]
# Note (!) the input data from treefam in REC and RME differs (see io tests)

# inf - input file from the selection above
# model - LCA, GMS, PG
# res - published results


@pytest.mark.parametrize("inf1, inf2, model, res",
                         [(inrecfiles[0], inrecfiles[1], 1, 7), (inrecfiles[0], inrecfiles[1], 2, 4),
                          (inrecfiles[0], inrecfiles[1], 3, 4), (inrecfiles[2], inrecfiles[3], 1, 47),
                          (inrecfiles[2], inrecfiles[3], 2, 46), (inrecfiles[2], inrecfiles[3], 3, 46)])
def test_rec(inf1, inf2, model, res):
    st = None
    gtrees = []
    gtrees.extend(readgtreefile(inf1))
    listtree = readgtreefile(inf2)
    if listtree:
        st = listtree[0]

    if model:
        for gt in gtrees:
            gt.set_lca_mapping(st)

    # generates intervals according to chosen model
    for gt in gtrees:
        if model == MOD_PASZEKGORECKI:
            genPaszekGoreckiIntervals(gt, st)
        elif model == MOD_GUIGO:
            genGMSIntervals(gt, st)
        elif model == MOD_LCA:
            genLCAIntervals(gt, st)

    out = rec(gtrees, st)
    assert out == res


@pytest.mark.parametrize("inf, model, res",
                         [(infiles[2], 1, 46), (infiles[2], 2, 45), (infiles[2], 3, 45),
                          (infiles[1], 1, 6), (infiles[1], 2, 4), (infiles[1], 3, 4),
                          (infiles[0], 1, 7), (infiles[0], 2, 4), (infiles[0], 3, 4)])
def test_rec(inf, model, res):
    gtrees, st = readintervalfile(inf)

    if model:
        for gt in gtrees:
            gt.set_lca_mapping(st)

    # generates intervals according to chosen model
    for gt in gtrees:
        if model == MOD_PASZEKGORECKI:
            genPaszekGoreckiIntervals(gt, st)
        elif model == MOD_GUIGO:
            genGMSIntervals(gt, st)
        elif model == MOD_LCA:
            genLCAIntervals(gt, st)

    out = rec(gtrees, st)
    assert out == res


@pytest.mark.parametrize("inf, res", [(infiles[0], 1), (infiles[1], 1), (infiles[2], 1)])
def test_rec_fhs(inf, res):
    gtrees, st = readintervalfile(inf)
    for gt in gtrees:
        gt.set_lca_mapping(st)
    for gt in gtrees:
        genFHSIntervals(gt, st)
    out = rec(gtrees, st)
    assert out == res
