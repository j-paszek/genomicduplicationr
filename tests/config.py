import os

cur_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    'data',
    )

ALL_INFILES = [
    # Following 3 files are from attached RME publication materials - more description see test_rme.py
    os.path.join(cur_dir, "RME", "in.txt"),                     # Guigo dataset (st)
    os.path.join(cur_dir, "RME", "in1.txt"),                    # Guigo dataset (st1)
    os.path.join(cur_dir, "RME", "in2.txt"),                    # Treefam dataset
    # Following files are from attaches UEC publication materials - more description see test_rec.py
    #       Note: publication focuses on input unrooted gene trees
    #       Guigo and Treefam datasets original files consists of rooted trees
    #       below are files from publication materials
    os.path.join(cur_dir, "REC", "guigo", "g.txt"),
    os.path.join(cur_dir, "REC", "guigo", "s1.txt"),
    os.path.join(cur_dir, "REC", "treefam", "g.txt"),
    os.path.join(cur_dir, "REC", "treefam", "s.txt"),
    #       however, during computation process trees were re-rooted
    #       according to the rules described in the publication
    #       In the publication are described results of EC score computed on altered files.
    #       Below files are generated by 'algorithm2.sh'
    #       (script version and 'urec' version - from sources from UEC publication materials)
    os.path.join(cur_dir, "REC", "guigo", "grootingssamplebeta.txt"),
    os.path.join(cur_dir, "REC", "treefam", "grootingssamplebeta.txt"),
    ]

MOD_LCA = 1  # LCA model
MOD_GUIGO = 2  # GMS model Guigo et al.
MOD_PASZEKGORECKI = 3  # PG model
MOD_FELLOWS = 4  # FHS model Fellows et al.
