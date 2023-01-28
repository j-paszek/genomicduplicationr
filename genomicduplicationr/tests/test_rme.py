import pytest
from genomicduplicationr.iomod import readintervalfile
from genomicduplicationr.rme import genFHSIntervals, genLCAIntervals, genGMSIntervals, genPaszekGoreckiIntervals, rme
from genomicduplicationr.rme_fhs import merfellows

MOD_LCA = 1  # LCA model
MOD_GUIGO = 2  # GMS model Guigo et al.
MOD_PASZEKGORECKI = 3  # PG model
MOD_FELLOWS = 4  # FHS model Fellows et al.

# Test corresponds to published results in
# Paszek, J., Gorecki, P., 2018. Efficient algorithms for genomic duplication models.
# IEEE/ACM Trans. Comput. Biol. Bioinform. 15 (5), 1515-1524.
# The results are presented in Figure 7 and in Figure 9
infiles = ["data/RME/in.txt", "data/RME/in1.txt", "data/RME/in2.txt"]


# inf - input file from the selection above
# model - LCA, GMS, PG
# res - published results


@pytest.mark.parametrize("inf, model, res", [(infiles[0], 1, 12), (infiles[0], 2, 7), (infiles[0], 3, 6),
                                             (infiles[1], 1, 9), (infiles[1], 2, 5), (infiles[1], 3, 5),
                                             (infiles[2], 1, 249), (infiles[2], 2, 243), (infiles[2], 3, 230)])
def test_rme(inf, model, res):
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

    out = rme(gtrees, st)
    assert out == res


@pytest.mark.parametrize("inf, res", [(infiles[0], 4), (infiles[1], 4)])
def test_rme_fhs(inf, res):
    gtrees, st = readintervalfile(inf)
    for gt in gtrees:
        gt.set_lca_mapping(st)
    for gt in gtrees:
        genFHSIntervals(gt, st)
    out = merfellows(gtrees, st)
    assert out == res