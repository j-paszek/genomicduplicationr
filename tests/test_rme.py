import pytest
from gdscore.iomod import readintervalfile
from gdscore.rme import genFHSIntervals, genLCAIntervals, genGMSIntervals, genPaszekGoreckiIntervals, rme
from gdscore.rme_fhs import merfellows
from tests.config import ALL_INFILES, MOD_LCA, MOD_GUIGO, MOD_PASZEKGORECKI

# All tests correspond to results published in
# Paszek, J., Gorecki, P., 2018. Efficient algorithms for genomic duplication models.
# IEEE/ACM Trans. Comput. Biol. Bioinform. 15 (5), 1515-1524.
# The results are presented in Figure 7 and in Figure 9
# Guigo dataset (st)    - ALL_INFILES[0] - REC score -  12,   7,   6,   4  for LCA, GMS, PG, FHS models
# Guigo dataset (st1)   - ALL_INFILES[1] - REC score -   9,   5,   5,   4  for LCA, GMS, PG, FHS models
# Treefam dataset       - ALL_INFILES[2] - REC score - 249, 243, 230       for LCA, GMS, PG      models

# inf - input files (see also description in config.py)
# model - LCA, GMS, PG
# res - published results
@pytest.mark.parametrize("inf, model, res",
                         [(ALL_INFILES[0], 1, 12), (ALL_INFILES[0], 2, 7), (ALL_INFILES[0], 3, 6),
                          (ALL_INFILES[1], 1, 9), (ALL_INFILES[1], 2, 5), (ALL_INFILES[1], 3, 5),
                          (ALL_INFILES[2], 1, 249), (ALL_INFILES[2], 2, 243), (ALL_INFILES[2], 3, 230)])
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


# inf - input files (see also description in config.py)
# model - FHS
# res - published results
@pytest.mark.parametrize("inf, res", [(ALL_INFILES[0], 4), (ALL_INFILES[1], 4)])
def test_rme_fhs(inf, res):
    gtrees, st = readintervalfile(inf)
    for gt in gtrees:
        gt.set_lca_mapping(st)
    for gt in gtrees:
        genFHSIntervals(gt, st)
    out = merfellows(gtrees, st)
    assert out == res
