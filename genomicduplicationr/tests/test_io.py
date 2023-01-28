import pytest
from genomicduplicationr import iomod

# Test works on input files stored in data/ that consist of gene trees and a species tree
# The files are loades, trees are generated and first gene tree and a species tree
# are converted back to string and compared with a corresponding part of an input file
infiles = ["data/RME/in.txt", "data/RME/in1.txt", "data/RME/in2.txt"]
inrecfiles = ["data/REC/guigo/g.txt", "data/REC/guigo/s1.txt",
              "data/REC/treefam/g.txt", "data/REC/treefam/s.txt"]
# ! NOTE ! : input gene trees mismatch for treefam in REC and RME
# the treefam from RME (data/RME/in2.txt)           - has 1274 trees; 1205 unique;
# the treefam from REC (data/REC/treefam/g.txt)     - has 1274 trees; 1206 unique; (838 trees match)


st1 = "(prot,(fung,((chlo,embr),(arth,((acoe,anne),(echi,(chon,(oste,(amph,(moll,((mamm,(aves,rept)),agna)))))))))))"
st2 = "(fung,(prot,((chlo,embr),((acoe,(arth,(anne,moll))),(echi,((agna,(chon,oste)),(amph,(mamm,(aves,rept)))))))))"
gt3 = "(((ORYSA,(ORYSA,(((ARATH,ARATH),ARATH),(ARATH,ARATH)))),(CIOIN,CAEEL)),((((ORYSA,ORYSA),((((ARATH,ARATH)," \
      "ARATH),ARATH),ARATH)),(ARATH,(((ARATH,ARATH),ARATH),(ARATH,ARATH)))),((((SCHPO,SCHPO),((YEAST,YEAST),(YEAST," \
      "YEAST))),((SCHPO,SCHPO),(YEAST,YEAST))),((((((((((((MOUSE,MOUSE),RAT),(((HUMAN,PANTR),MACMU),((BOVIN,BOVIN)," \
      "((CANFA,CANFA),CANFA)))),MONDO),XENTR),XENTR),(BRARE,((BRARE,BRARE),(TETNG,GASAC)))),((((((((HUMAN,PANTR)," \
      "MACMU),(BOVIN,CANFA)),(MOUSE,(RAT,RAT))),MONDO),CHICK),(XENTR,XENTR)),(BRARE,(TETNG,(GASAC,GASAC))))),(DROME," \
      "(AEDAE,(ANOGA,ANOGA)))),(((CAEEL,CAEEL),CAEEL),CAEBR)),SCHMA),(((((((((HUMAN,MACMU),CANFA),(MOUSE,RAT))," \
      "MONDO),XENTR),(BRARE,TETNG)),CIOIN),((DROME,(AEDAE,ANOGA)),(CAEEL,CAEBR))),SCHMA)))))"
st3 = "((((((((AEDAE,ANOGA),(DROME,DROPS)),APIME),SCHMA),(((((((BOVIN,CANFA),(((HUMAN,PANTR),MACMU),(MOUSE,RAT)))," \
      "MONDO),CHICK),XENTR),(BRARE,((TETNG,FUGRU),(ORYLA,GASAC)))),CIOIN)),(CAEBR,CAEEL)),(YEAST,SCHPO)),(ARATH,ORYSA))"


@pytest.mark.parametrize("infile,tgt,tst", [(infiles[0], "(((aves,mamm),arth),prot)", st1),
                                            (infiles[1], "(((aves,mamm),arth),prot)", st2),
                                            (infiles[2], gt3, st3)])
def test_readintervalfile(infile, tgt, tst):
    gtrees, st = iomod.readintervalfile(infile)
    assert str(gtrees[0]) == tgt
    assert str(st) == tst


@pytest.mark.parametrize("inf_all,inf_gt,inf_st", [(infiles[0], inrecfiles[0], inrecfiles[1])])
def test_variuos_inputs(inf_all, inf_gt, inf_st):
    gtrees, st = iomod.readintervalfile(inf_all)
    gts = iomod.readgtreefile(inf_gt)
    st1 = iomod.readgtreefile(inf_st)

    for i in range(len(gtrees)):
        assert str(gtrees[i]) == str(gts[i])
    assert str(st) == str(st1[0])


@pytest.mark.parametrize("inf_all,inf_gt,inf_st", [(infiles[2], inrecfiles[2], inrecfiles[3])])
def test_variuos_inputs_tf(inf_all, inf_gt, inf_st):
    gtrees, st = iomod.readintervalfile(inf_all)
    gts = iomod.readgtreefile(inf_gt)
    st1 = iomod.readgtreefile(inf_st)

    assert len(gtrees) == len(gts)
    assert str(st) == str(st1[0])
