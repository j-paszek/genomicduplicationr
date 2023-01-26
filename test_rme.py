import pytest
import rme

ModLCA=1                # LCA model
ModGuigo=2              # GMS model Guigo et al.
ModPaszekGorecki=3      # PG model
ModFellows=4            # FHS model Fellows et al.

# FIXME add tests for Fellows model
# Test corresponds to published results in
# Paszek, J., Gorecki, P., 2018. Efficient algorithms for genomic duplication models.
# IEEE/ACM Trans. Comput. Biol. Bioinform. 15 (5), 1515-1524.
# The results are presented in Figure 7 and in Figure 9
infiles = ["data/in.txt","data/in1.txt","data/in2.txt"]
# inf - input file from the selection above
# model - LCA, GMS, PG
# res - published results
@pytest.mark.parametrize("inf,model,res", [(infiles[0],1,12),(infiles[0],2,7),(infiles[0],3,6),
                                           (infiles[1],1,9),(infiles[1],2,5),(infiles[1],3,5),
                                           (infiles[2],1,249),(infiles[2],2,243),(infiles[2],3,230)])
def test_mer(inf,model,res):
    gtrees, st = rme.readintervalfile(inf)

    if model:
        for gt in gtrees:
            gt.setlcamapping(st)

    for gt in gtrees:
        if model==ModPaszekGorecki:
            rme.genPaszeGoreckiIntervals(gt,st)
        elif model==ModGuigo:
            rme.genGMSIntervals(gt,st)
        elif model==ModLCA:
            rme.genLCAIntervals(gt,st)

    out = rme.mer(gtrees,st)
    assert out == res



