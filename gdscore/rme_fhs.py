import sys
from itertools import combinations, chain
from gdscore.rme import expandintervals, rme, meropt, get_left, get_right, is_dup


def speciestreescore(n):
    if n.leaf():
        s = n.clusterleaf
    else:
        s = "(" + speciestreescore(n.r) + "," + speciestreescore(get_left(n)) + ")"
    if n.scorek > 0:
        s = s + " k=%d " % n.scorek
    return s


def ppscores(st):
    total = sum(n.scorek for n in st.nodes)
    return speciestreescore(st.root)+" treename='ME=%d'" % total


def all_subsets(ss):
    return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))


def merfellows(gtrees, st, verbose="", printstreewithscores=0):
    lastopt = ''
    if verbose.find("4") != -1:
        print("@ Debug of Fellows model")
    specnodes = [n for g in gtrees for n in g.root.nodes()
                 if not n.leaf() and (n.lcamap != get_left(n).lcamap and n.lcamap != get_right(n).lcamap)]
    if verbose.find("4") != -1:
        print("@ Number all speciation nodes before preprocessing", len(specnodes))

    r = []
    for s in specnodes:
        hasdup = 0
        for n in s.nodes():
            if not n.leaf() and (is_dup(n)):
                hasdup = 1
                break
        if hasdup:
            r.append(s)
    specnodes = r

    if verbose.find("4") != -1:
        print("@ List of all speciation nodes", specnodes)
        print("@ Number of all speciation nodes", len(specnodes))

    # the following commented legacy code is related to the treefam dataset
    # some nodes were analysed (see the publication Paszek, Gorecki, 2018.
    # Efficient algorithms for genomic duplication models TCBB)
    # and manually excluded to make the computation feasible
    autodups = []
    # autodupsnum=[66, 85, 91, 106, 186, 208, 204, 69, 243, 247, 254, 258, 253, 270, 268, 307, 293, 300]
    # # HACK dla 490
    # for i,g in enumerate(gtrees):
    #     for n in g.nodes:
    #         if n.num in autodupsnum: 
    #              subnodes.extend(n.nodes())
    # autodups=[ n for n in specnodes if n in subnodes ]
    # specnodes=[ n for n in specnodes if n not in subnodes ]
    # #HACK dla 1257
    # for i,g in enumerate(gtrees):
    #     for n in g.nodes:
    #         if n.lcamap.num==39: # should be <= 39
    #             subnodes.extend(n.nodes())
    # autodups=[ n for n in specnodes if n in subnodes ]
    # specnodes=[ n for n in specnodes if n not in subnodes ]
    # #END OF HACK

    if verbose.find("2") != -1:
        print("Found ", len(specnodes), "speciation nodes")
        print("@ List of all speciation nodes", specnodes)

    maxdup = 0
    alldupscore = 0
    pompom = None
    pom = None

    for i, g in enumerate(gtrees):
        maxh = maxc = -1
        for n in g.root.nodes():
            if n.leaf():
                c = 0
        
                while n:
                    if not n.leaf() and (is_dup(n)):
                        c = c+1
                    n = n.parent
            if c > maxdup:
                pom = g
            maxdup = max(maxdup, c)
            maxc = max(c, maxc)
        if alldupscore < g.height():
            pompom = g
        alldupscore = max(alldupscore, g.height())
        if verbose.find("6") != -1:
            print(i, maxh, maxc)

    if verbose.find("2") != -1:
        print("Bottom limit for ME score:", maxdup)
        print("Bottom limit tree:", pom)
        print("AllDupScore:", alldupscore)
        print("AllDup tree:", pompom)

    if verbose.find("6") != -1:
        # mer.py -F -j 1257,1 -v6 -e in.txt
        # in treefam

        if len(gtrees) > 10:
            print("too many trees")
            sys.exit(-1)  # stop    

        def ptree(t):
            if t.leaf(): 
                t.lf = 1
                return
            if is_dup(t):
                t.lf = 0
                t.lcamap.ds |= 1
            elif t in specnodes:
                t.lcamap.ds |= 2
                t.lf = 0
            else:
                t.lf = 1
            if not t.lf:
                ptree(get_left(t))
                ptree(get_right(t))

        def ppniceold(t, isst):  # up to two marks in S
            l = ''
            dp = 0
            if not isst:
                if t.leaf(): 
                    return "lw=0.2 s=\"%s\" :0.0" % t.clusterleaf
                if is_dup(t):
                    l = "mark=d%d nn=\"%d\"" % (t.lcamap.num, t.num)
                    dp = 1
                if t in specnodes:
                    l = "mark=s%d nn=\"%d\"" % (t.lcamap.num, t.num)
                                
                if not l:
                    # return "mark=94 :0.3"#%(t.num,t.num)
                    return ":0.0"  # +("mark=100" if t.lf else "")
                    # 2023 - comment next 3 lines
                    # l="mark=s%d "%t.lcamap.num
                    # t.lcamap.ds|=2
                    # return l+"s=''" #\"%d\""%t.lcamap.num
            
                # x=" s=\"!R!%d \""%t.lcamap.num
                x = ''
            
                if t.lcamap.leaf() and dp == 1:
                    if get_right(t).stheight > get_left(t).stheight:
                        return "("+ppniceold(get_right(t), isst)+")"+l
                    return "("+ppniceold(get_left(l), isst)+")"+l
            else:                 
                if t.ds == 3:
                    x = " mark=[d%d,s%d] " % (t.num, t.num)
                elif t.ds == 2:
                    x = " mark=[s%d] " % t.num
                elif t.ds == 1:
                    x = " mark=[d%d] " % t.num
                else:
                    x = ''
                x = x+" nn=\"%d\"" % t.num
                
                if t.leaf():
                    return x+" s=\"%s\"" % t.clusterleaf
                # else: x=x+" s=\"!R!%d \""%(t.num)
            if get_left(t).lf:
                return "("+ppniceold(get_right(t), isst)+") "+l+x
            elif get_right(t).lf:
                return "("+ppniceold(get_left(t), isst)+") "+l+x
            return "("+ppniceold(get_left(t), isst)+","+ppniceold(get_right(t), isst)+") "+l+x

        def ppnice(t, isst):  # single mark in S
            l = ''
            dp = 0
            maxmark = 29*3
            if not isst:
                if t.leaf(): 
                    return "lw=0.2 s=\"%s\" :0.0" % t.clusterleaf
                if is_dup(t):
                    l = "mark=d%d nn=\"%d\"" % (t.lcamap.num, t.num)
                    t.lcamap.ds |= 1
                    dp = 1
                if t in specnodes:
                    l = "mark=s%d nn=\"%d\"" % (t.lcamap.num, t.num)
                    t.lcamap.ds |= 2
            
                if not l:
                    # return "mark=94 :0.3"#%(t.num,t.num)
                    return ":0.0"  # +("mark=100" if t.lf else "")
                    
                x = ''
                if t.lcamap.leaf() and dp == 1:
                    if get_right(t).stheight > get_left(t).stheight:
                        return "("+ppnice(get_right(t), isst)+")"+l
                    return "("+ppnice(get_left(t), isst)+")"+l
            else:                 
                if t.ds == 3:
                    x = " mark=[d%d,s%d] " % (t.num, t.num)
                elif t.ds == 2:
                    x = " mark=s%d " % t.num
                elif t.ds == 1:
                    x = " mark=d%d " % t.num
                else:
                    x = ''
                x = x+" nn=\"%d\"" % t.num
                
                if t.leaf():
                    return x+" s=\"%s\"" % t.clusterleaf
                # else: x=x+" s=\"!R!%d \""%(t.num)
            if get_left(t).lf:
                return "("+ppnice(get_right(t), isst)+") "+l+x
            elif get_right(t).lf:
                return "("+ppnice(get_left(t), isst)+") "+l+x
            return "("+ppnice(get_left(t), isst)+","+ppnice(get_right(t), isst)+") "+l+x
        f = open("dup.gse", "w")
        f.write("""
.run -dm
leafline=5
treestyle=2
lcamapping=0
evalpos=3
rotateleaves=90
gene.scale=(0.6,.6)
species.scale=(0.7,.7)
multimarkalign='c'
mapsep=-5
lw=0.7
gene.leaflabellen=1
species.leaflabellen=55
species.nodelabels=["s"]
gene.nodelabels=["z"]
markscheme="rgb W G"
showleaves=1
mapbymarks=1
marksize=1.0
""")
        for n in st.nodes: 
            n.ds = 0
            n.lf = 0
            n.sp = n.dp = 0
        for g in gtrees:
            for n in g.nodes:
                n.lf = 0
            ptree(g.root)
            cnt = 0
            for n in st.nodes:
                if n.ds or n.sp: 
                    n.marknum = cnt
                    cnt += 1
            f.write("&g %s\n" % ppnice(g.root, 0))
        f.write("&s %s\n" % ppnice(st.root, 1))
        cnt = 0
        maxmark = 29*3
        for n in st.nodes: 
            if n.ds & 1:
                f.write("d%d=%d\n" % (n.num, n.marknum))
            if n.ds & 2:
                f.write("s%d=%d\n" % (n.num, n.marknum+maxmark))

        # print "Dupnodesmap:",cnt,"Specnodesmap",scnt-23*3
        f.close()
        print("File dup.gse saved")
        quit()

    mescore = sum(len(gt.nodes) for gt in gtrees)+1
    first = 1
    cnt = 0
    for ss in all_subsets(specnodes):         # for ss in [ specnodes ]:
        if verbose.find("4") != -1:
            print("@ Subset of speciation nodes enabled for a change into duplication nodes", ss)
        ss = list(ss)
        s = [n for n in specnodes if n not in ss ]
        # Creating intervals
        for gt in gtrees:
            for n in gt.nodes:
                n.interval = None
                if n.leaf(): continue
                if is_dup(n):
                    n.interval = [n.lcamap, None]

        for n in s:
            n.interval = [n.lcamap, None]

        # HACK
        for n in autodups:
            n.interval = [n.lcamap, None]

        for gt in gtrees:
            expandintervals(gt, st)

        if verbose.find("4") != -1:
            print("@ List of all intervals (spec->dup node intervals are extended)")
            for gt in gtrees:            
                for n in gt.nodes:                
                    if n.interval:                    
                        print("Node ", n, " interval ", n.interval)
        
        if first:
            m = mescore = rme(gtrees, st, verbose)
            first = 0
            if verbose.find("1") != -1:
                print("Current score", m)
        else:
            m = meropt(gtrees, st, mescore, verbose)
        if m > -1:
            if verbose.find("4") != -1:
                print("%d. MEcurrent" % cnt, m, mescore)
            if mescore > m:
                if verbose.find("1") != -1:
                    oldverbose, verbose = verbose, ""
                    print("*"*80)
                    print(cnt, "REPEATED OPTIMAL COMP FOR ", m)
                meropt(gtrees, st, -1, verbose)
                lastopt = ppscores(st)
                if verbose.find("1") != -1:
                    print("*"*80)
                    verbose = oldverbose
                    print("Current score", m)
            mescore = min(mescore, m)
        cnt = cnt+1
        if verbose.find("1") != -1:
            if cnt % 1000 == 0:
                print(cnt, "variants processed; current min", mescore, "last score", m)
        # if cnt==5000: quit()
    if printstreewithscores:
        print("&s", lastopt)
    else:
        if verbose.find("1") != -1:
            print("MEscore", mescore)
    return mescore
