from treeop import Tree, str2tree
import itertools


def rme(gtrees, st, verbose=""):
    return meropt(gtrees, st, -1, verbose)


def meropt(gtrees, st, prevminscore, verbose):
    dup = []
    stnodespostorder = st.root.nodes()
    gtnodespostorder = list(itertools.chain.from_iterable(gt.root.nodes() for gt in gtrees))

    for g in gtnodespostorder:
        g.active = bool(g.interval)
        if g.active:
            dup.append(g)
        g.visited = None

    if verbose.find("1") != -1:
        print("Duplication Nodes", dup)
    # print stnodespostorder
    
    # set interval nodes in S 
    for s in stnodespostorder:
        s.topinterval = []
        s.botinterval = []
        s.bt = []
        s.duproots = []  # Gamma_t's
        s.scorek = 0  # clean k

    for g in dup:
        g.interval[0].botinterval.append(g)
        g.interval[1].topinterval.append(g)

    # Nodes of T (by topintervals)
    tdup = [s for s in stnodespostorder if s.topinterval]

    # leaves of duplication forests
    ldup = [g for g in gtnodespostorder if g.active and (g.leaf() or (not g.l.active and not g.r.active))]

    # Set B sets
    def setnearestt(n, t):
        if n.topinterval:
            n.parentt = t
            t = n
        n.nearestt = t
        if n.leaf():
            return
        setnearestt(n.l, t)
        setnearestt(n.r, t)

    setnearestt(st.root, None)

    for l in ldup:
        l.interval[0].nearestt.bt.append(l)

    # find Gamma's
    for d in dup:
        # not rootT || parent nie jest dupl || inny max dup parenta )
        if not d.parent or not d.parent.active or (d.parent.active and d.parent.interval[1] != d.interval[1]):
            d.interval[1].duproots.append(d)   

    if verbose.find("1") != -1:
        print("Nodes of T")
        for t in tdup: 
            print(t)
            print("   Bt=", t.bt, [n.num for n in t.bt])
            print("   Gammat=", t.duproots)

    def travduproots(r):
        if r.leaf():
            return 1
        
        h1 = h2 = 0
        if r.l.active:
            h1 = travduproots(r.l)
        if r.r.active:
            h2 = travduproots(r.r)
        return 1+max(h1, h2)

    if verbose.find("1") != -1:
        print("Dup Leaves", ldup)

    mescore = 0

    # Main loop
    for t in tdup:
        if verbose == 1:        
            print("="*80)
            print("Processing", t)
            print("   Bt=", t.bt)
            print("   Gammat=", t.duproots)

        k = -1
        for r in t.duproots:
            if r.active:
                k = max(k, travduproots(r))
        if verbose.find("1") != -1:
            print("   k=", k)

        if k > 0:
            mescore += k
            t.scorek = k
        else:
            t.scorek = 0

        cand = []
        if k < 0:
            if verbose.find("1") != -1:
                for l in t.bt:
                    print("Moving from orhpan dup leaf", l, l.active)
            t.parentt.bt.extend(t.bt)

        else:
            # processing leaves from Bt, k>0 case
            for l in t.bt:
                if verbose.find("1") != -1:
                    print("Removing from dup leaf", l.num, l, l.active)
                curk = 1
                newcand = None
                while curk <= k and not newcand:
                    l.visited = t
                    if verbose.find("1") != -1:
                        print("  Marking", l.num, l)
                    if not l.parent:
                        break  # gene tree root

                    l.h = curk
                    sib = l.sibling()
                    l = l.parent
                    if verbose.find("1") != -1:
                        print("  checking", l, l.num if l else "")  # l.interval
                    if not l.active:
                        if verbose.find("1") != -1:
                            print("  not active")
                        break
                    if verbose == 1:
                        print(l.interval)
                    
                    if l.interval[0].nearestt.lca(t) != t or curk == k:
                        newcand = l
                        if verbose.find("1") != -1:
                            print(" cand. found", l)
                        # stop
                    elif not sib.active:
                        curk += 1
                        if verbose.find("1") != -1:
                            print("  sib is not active - cont", sib)
                    elif sib.visited == t:
                        curk = max(curk, sib.h)+1
                    else:
                        if verbose.find("1") != -1:
                            print("  sib not visited yet", sib)
                        break

                if newcand:
                    cand.append(newcand)

            # Clean active nodes
            for l in t.bt:
                while l and l.visited == t:
                    l.active = 0
                    l = l.parent

        if verbose.find("1") != -1:
            print("Candidates", cand)

        # processing bt candidates
        for c in cand:
            if not c.l.active and not c.r.active:
                if c.interval[0].nearestt.lca(t) == t:
                    if verbose.find("1") != -1:
                        print(" cand ", c, "moved to (parentt)", t.parentt)
                    t.parentt.bt.append(c)
                    # ??? Nono (root T)?
                else:
                    # if k<0: print "Wrr"
                    c.interval[0].nearestt.bt.append(c)
                    if verbose.find("1") != -1:
                        print(" cand ", c, "moved to (>parentt)", c.interval[0].nearestt)
                             
        if verbose.find("1") != -1:
            print("Remaining dupl:")
        for d in dup:
            if d.active:
                if verbose.find("1") != -1:
                    print(" ", d.num, d)  # ,hasattr(d,'visited')
        if prevminscore >= 0 and mescore >= prevminscore:
            if verbose.find("4") != -1:
                print("Skipping score ", mescore, "vs", prevminscore)
                return -1
            return mescore 
    if verbose.find("2") != -1:
        print("MEscore", mescore)
    return mescore
            

def genLCAIntervals(gt, st):
    for n in gt.nodes:
        n.interval = None
        if n.leaf():
            continue
        if n.lcamap == n.l.lcamap or n.lcamap == n.r.lcamap:
            n.interval = [n.lcamap, n.lcamap]


def genPaszekGoreckiIntervals(gt, st):
    for n in gt.nodes:
        n.interval = None
        if n.leaf():
            continue
        if n.lcamap == n.l.lcamap or n.lcamap == n.r.lcamap:
            n.interval = [n.lcamap, None]

    expandintervals(gt, st)


def expandintervals(gt, st):
    for n in gt.nodes:
        if n.interval:        
            c = n.parent
            while c and c.interval:
                c = c.parent
            if not c: 
                n.interval[1] = st.root
                continue
            top = n.interval[0]
            while top.parent != c.lcamap:
                top = top.parent
            n.interval[1] = top


def genGMSIntervals(gt, st):
    for n in gt.nodes:
        n.interval = None
        if n.leaf():
            continue
        if n.lcamap == n.l.lcamap or n.lcamap == n.r.lcamap:
            n.interval = [n.lcamap, None]

    for n in gt.nodes:
        if not n.interval:
            continue
        # if n.num==19: print "ASDADA",n
        if not n.parent:  # root
            n.interval[1] = st.root
            continue        

        if n.parent.interval and n.lcamap == n.parent.lcamap:
            n.interval[1] = n.interval[0]
            continue
            
        top = n.lcamap

        while top.parent != n.parent.lcamap:
            # if n.num==19:
            #    print "-----",top,n.lcamap,n.parent.lcamap
            top = top.parent
        # if n.num==19:
            # print "--LAST",top,n.lcamap,n.parent.lcamap
        n.interval[1] = top


def ppgsevoloutput(gtrees, st):
    f = open("mer.gse", "w")
    f.write("#!/home/gorecki/Dropbox/projects/htree/gsevol.py -x\n")
    for i, gt in enumerate(gtrees):
        f.write("""
    
####################################################
# MAPPINGS FOR DL

#exec:Mappings for DL: -dDPm
outputfile='%d.map'

#colonsep

mapprops=[
('green',0,0,0.5,'') ,
('blue',0,0,0.5,'') ,
('red',0,1,0.5,'') ,
('green',0,0,0.5,'') ,
('blue',0,0,0.5,'') ,
];

#colonsep

lw=0.7
resolution=2
crval=0.45
dotsep=1.4
crvbe=0.15
arrowlength=0.2
mapleaves=0
mapinternal=1
rootedgelen=0.3
gtscaleto=(30,100)
stscaleto=(50,10)
scale=(0.6,0.7)
mappict=1
mapduponly=1
sep=40
margin=15
maplw=0.4
leaflabels=1
nodelabels=['total','leaflabel']
extsyntax=1
rotateleaves=0
fillstyle=SOLID
mapdot=1
treestyle=1
evalpos=4

tableonly=1
table=[20,4] 

&g %s
&s %s 


""" % (i, gt, st))

    f.close()

    
def readgtreefile(gtreefile):
    gtrees = []
    for g in open(gtreefile).readlines():
        if g and g[0] != "#":
            gtrees.append(Tree(str2tree(g.strip())))

    return gtrees
