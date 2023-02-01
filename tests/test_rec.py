import pytest
from genomicduplicationr.iomod import readgtreefile, readintervalfile
from genomicduplicationr.rme import genFHSIntervals, genLCAIntervals, genGMSIntervals, genPaszekGoreckiIntervals
from genomicduplicationr.rec import rec
from tests.config import ALL_INFILES, MOD_LCA, MOD_GUIGO, MOD_PASZEKGORECKI

# Test 1-4 corresponds to published results in
# J. Paszek and P. Gorecki, “Genomic duplication problems for unrooted gene trees,”
# BMC Genomics, vol. 17, no. 1, pp. 165–175, 2016.
#           The results correspond to values presented in table 1
# Note (!)  The paper focus on unrooted variant of the problem, hence,
#           the input gene trees reported in table 1 were rerooted. The details are described in the paper.
#           To recreate the rootings use "algorithm2.sh" (see https://www.mimuw.edu.pl/~jpaszek/uec.php) or use
#           attached files (in 'grootingssamplebeta.txt' are rerooted trees from corresponding file 'g.txt.)
# Note (!)  the input data from treefam in REC and RME differs (see io tests)
# Guigo dataset     - ALL_INFILES[7] (gtrees) ALL_INFILES[4] (st) - REC score -   5,  4  for GMS, PG  models
# Treefam dataset   - ALL_INFILES[8] (gtrees) ALL_INFILES[6] (st) - REC score -  45, 45  for GMS, PG  models


# inf - input files (see also description in config.py)
# model - LCA, GMS, PG
# res - results - first 4 tests correspond to published results; next - other tests based on not rerooted gene trees
@pytest.mark.parametrize("inf1, inf2, model, res",
                         [(ALL_INFILES[7], ALL_INFILES[4], 2, 5), (ALL_INFILES[7], ALL_INFILES[4], 3, 4),
                          (ALL_INFILES[8], ALL_INFILES[6], 2, 45), (ALL_INFILES[8], ALL_INFILES[6], 3, 45),
                          (ALL_INFILES[3], ALL_INFILES[4], 1, 7), (ALL_INFILES[3], ALL_INFILES[4], 2, 4),
                          (ALL_INFILES[3], ALL_INFILES[4], 3, 4), (ALL_INFILES[5], ALL_INFILES[6], 1, 47),
                          (ALL_INFILES[5], ALL_INFILES[6], 2, 46), (ALL_INFILES[5], ALL_INFILES[6], 3, 46)])
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


# inf - input files (see also description in config.py)
# model - LCA, GMS, PG
# res - results for RME files
@pytest.mark.parametrize("inf, model, res",
                         [(ALL_INFILES[2], 1, 46), (ALL_INFILES[2], 2, 45), (ALL_INFILES[2], 3, 45),
                          (ALL_INFILES[1], 1, 6), (ALL_INFILES[1], 2, 4), (ALL_INFILES[1], 3, 4),
                          (ALL_INFILES[0], 1, 7), (ALL_INFILES[0], 2, 4), (ALL_INFILES[0], 3, 4)])
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


# inf - input files (see also description in config.py)
# model - FHS
# res - for unrestricted model all duplications should cluster in one
@pytest.mark.parametrize("inf, res", [(ALL_INFILES[0], 1), (ALL_INFILES[1], 1), (ALL_INFILES[2], 1)])
def test_rec_fhs(inf, res):
    gtrees, st = readintervalfile(inf)
    for gt in gtrees:
        gt.set_lca_mapping(st)
    for gt in gtrees:
        genFHSIntervals(gt, st)
    out = rec(gtrees, st)
    assert out == res
