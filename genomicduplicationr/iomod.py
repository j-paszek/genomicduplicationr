from genomicduplicationr.treeop import Tree, str2tree


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


def savegsi(gtrees, st, outputfile):
    f = open(outputfile, "w")
    f.write("%d\n" % len(gtrees))
    for gt in gtrees:
        f.write(str(gt) + "\n")

    f.write("\n#Species tree\n" + str(st) + "\n\n")

    for n in st.root.nodes():
        f.write("#%d %s\n" % (n.num, n))

    for i, gt in enumerate(gtrees):
        f.write("\n")
        f.write("#"*10)
        f.write(" Tree nr %d\n" % i)
        for g in gt.root.nodes():
            if g.interval:
                f.write("#%d %s\n" % (g.num, g))
                f.write("%d;%d;%d;%d\n" % (i, g.num, g.interval[0].num, g.interval[1].num))
    f.close()
    print("File %s saved" % outputfile)


def readgtreefile(gtreefile):
    gtrees = []
    for g in open(gtreefile).readlines():
        if g and g[0] != "#":
            gtrees.append(Tree(str2tree(g.strip())))

    return gtrees
