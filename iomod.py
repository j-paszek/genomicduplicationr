from treeop import Tree, str2tree

verbose = ""

def readintervalfile(filename):
    t = [l.strip() for l in open(filename, "r").readlines() if len(l.strip()) and l.strip()[0] != '#']

    if t[0][0].isdigit():
        gtnum = int(t[0])
        offset = 1
    else:
        gtnum = 1
        offset = 0

    gtrees = []

    for i in range(gtnum):
        if verbose == 3:
            print("Processing line ", t[i + offset])
        gtrees.append(Tree(str2tree(t[i + offset])))
    st = Tree(str2tree(t[i + 1 + offset]))

    stnodespostorder = st.root.nodes()

    oldgtnum = -1

    for i in range(i + 2 + offset, len(t)):

        gtnum, gc, sc, l = t[i].replace("+", " +").split(";")
        l = int(l)
        gtnum = int(gtnum)

        if oldgtnum != gtnum:
            gt = gtrees[gtnum]
            gtnodespostorder = gt.root.nodes()
            oldgtnum = gtnum

        if gc[0].isdigit():
            g = gtnodespostorder[int(gc)]
            s = stnodespostorder[int(sc)]
            l = stnodespostorder[l]

        else:
            gc = gc.strip().split()
            sc = sc.strip().split()

            g = gt.findnodeplus(gc)
            s = st.findnodeplus(sc)

        if not g:
            raise Exception("Incorrect interval cluster " + str(gc))
        if not s:
            raise Exception("Incorrect interval cluster " + str(sc))

        g.setinterval(s, l)

    return gtrees, st
