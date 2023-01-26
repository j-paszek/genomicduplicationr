#!/usr/bin/python

import sys
from treeop import Tree, str2tree, node2label
import getopt,itertools

verbose=""
printstreewithscores=0

def readintervalfile(filename):

    try:
        t=[ l.strip() for l in open(filename,"r").readlines() if len(l.strip()) and l.strip()[0]!='#' ]
    except IOError as e:
        print filename," I/O error({0}): {1}".format(e.errno, e.strerror)
        quit()
    
    if t[0][0].isdigit():
        gtnum=int(t[0])
        offset=1
    else: 
        gtnum=1
        offset=0
    
    gtrees=[]

    for i in xrange(gtnum):
        if verbose.find("0") <> -1:
            print "Processing line ", t[i+offset]
        gtrees.append(Tree(str2tree(t[i+offset])))
    st=Tree(str2tree(t[i+1+offset]))

    stnodespostorder=st.root.nodes()

    oldgtnum=-1
   
    for i in xrange(i+2+offset,len(t)):
    
        gtnum,gc,sc,l=t[i].replace("+"," +").split(";")
        l=int(l)
        gtnum=int(gtnum)

        if oldgtnum!=gtnum:
            gt=gtrees[gtnum]
            gtnodespostorder=gt.root.nodes()
            oldgtnum=gtnum
        
        if gc[0].isdigit():
            g=gtnodespostorder[int(gc)]
            s=stnodespostorder[int(sc)]
            l=stnodespostorder[l]

        else:
            gc=gc.strip().split()
            sc=sc.strip().split()


            g=gt.findnodeplus(gc)
            s=st.findnodeplus(sc)


        if not g:
            raise Exception("Incorrect interval cluster "+str(gc))
        if not s:
            raise Exception("Incorrect interval cluster "+str(sc))

        g.setinterval(s,l)

        
    return gtrees,st


def readgs(filename):

    t=[ l.strip() for l in open(filename,"r").readlines() if len(l.strip()) and l.strip()[0]!='#' ]
    if verbose.find("0") <> -1:
        print t
    #print str2tree(t[0])
    return Tree(str2tree(t[0])),Tree(str2tree(t[1]))

def savegsi(gtrees,st,outputfile):
    f=open(outputfile,"w")
    f.write("%d\n"%len(gtrees))
    for gt in gtrees: 
        f.write(str(gt)+"\n")

    f.write("\n#Species tree\n"+str(st)+"\n\n")

    for n in st.root.nodes():
        f.write("#%d %s\n"%(n.num,n))

    for i,gt in enumerate(gtrees):         
        f.write("\n")
        f.write("#"*10)
        f.write(" Tree nr %d\n"%i)
        for g in gt.root.nodes():
            if g.interval:
                f.write("#%d %s\n"%(g.num,g))
                f.write("%d;%d;%d;%d\n"%(i,g.num,g.interval[0].num,g.interval[1].num))


    f.close()
    print "File %s saved"%outputfile


def mer(gtrees,st):
    return meropt(gtrees,st,-1)

def meropt(gtrees,st,prevminscore):

    dup=[]
    
    stnodespostorder=st.root.nodes()
    gtnodespostorder=list(itertools.chain.from_iterable(gt.root.nodes() for gt in gtrees))

    for g in gtnodespostorder:
        g.active=bool(g.interval)    
        if g.active: dup.append(g)
        g.visited=None

    if verbose.find("2") <> -1:    
        print "Duplication Nodes",dup

    #print stnodespostorder
    
    # set interval nodes in S 
    for s in stnodespostorder:
        s.topinterval=[]
        s.botinterval=[]
        s.bt=[]
        s.duproots=[]  # Gamma_t's
        s.scorek=0  # clean k

    for g in dup:
        g.interval[0].botinterval.append(g)
        g.interval[1].topinterval.append(g)

    # Nodes of T (by topintervals)
    tdup=[ s for s in stnodespostorder if s.topinterval]

    # leaves of duplication forests
    ldup=[ g for g in gtnodespostorder if g.active and (g.leaf() or (not g.l.active and not g.r.active)) ]


    # Set B sets
    def setnearestt(n,t):
        if n.topinterval:
            n.parentt=t 
            t=n
        n.nearestt=t
        if n.leaf(): return
        setnearestt(n.l,t)
        setnearestt(n.r,t)

    setnearestt(st.root,None)

    for l in ldup:
        l.interval[0].nearestt.bt.append(l)

    # find Gamma's
    for d in dup:
        # not rootT || parent nie jest dupl || inny max dup parenta )
        if not d.parent or not d.parent.active or (d.parent.active and d.parent.interval[1]!=d.interval[1]):
            d.interval[1].duproots.append(d)   


    if verbose.find("2") <> -1:
        print "Nodes of T"
        for t in tdup: 
            print t
            print "   Bt=",t.bt
            print "   Gammat=",t.duproots

    def travduproots(r):
        if r.leaf(): return 1
        
        h1=h2=0
        if r.l.active: h1=travduproots(r.l)
        if r.r.active: h2=travduproots(r.r)
        return 1+max(h1,h2)

    if verbose.find("2") <> -1:    
        print "Dup Leaves",ldup

    mescore=0


    # Main loop
    for t in tdup:
        if verbose.find("2") <> -1:        
            print "="*80        
            print "Processing",t        
            print "   Bt=",t.bt        
            print "   Gammat=",t.duproots

        
        k=-1        
        for r in t.duproots:
            if r.active:
                k=max(k,travduproots(r))
        if verbose.find("2") <> -1:        
            print "   k=",k

        if k>0: 
            mescore+=k
            t.scorek=k
        else: t.scorek=0

    
        cand=[]
        if k<0:
            if verbose.find("2") <> -1:            
                for l in t.bt:
                    print "Moving from orhpan dup leaf",l,l.active
            
            t.parentt.bt.extend(t.bt)

        else:
            #processing leaves from Bt, k>0 case 
            
            for l in t.bt:
                if verbose.find("2") <> -1:                 
                    print "Removing from dup leaf",l.num,l,l.active
                curk=1
                newcand=None
                while curk<=k and not newcand:
                    l.visited=t
                    if verbose.find("2") <> -1:                    
                        print "  Marking",l.num,l
                    if not l.parent: break # gene tree root 

                    l.h=curk
                    sib=l.sibling()
                    l=l.parent
                    if verbose.find("2") <> -1:
                        print "  checking",l,l.num if l else "", # l.interval
                    if not l.active:
                        if verbose.find("2") <> -1: 
                            print "  not active"
                        break
                    if verbose == 1:
                        print l.interval
                    
                    if l.interval[0].nearestt.lca(t)!=t or curk==k:
                        newcand=l
                        if verbose.find("2") <> -1: print " cand. found",l
                        # stop
                    elif not sib.active:
                        curk+=1
                        if verbose.find("2") <> -1: print "  sib is not active - cont",sib
                    elif sib.visited==t:
                        curk=max(curk,sib.h)+1                
                    else:
                        if verbose.find("2") <> -1:                        
                            print "  sib not visited yet",sib
                        break

                if newcand: cand.append(newcand)

            # Clean active nodes
            for l in t.bt:
                while l and l.visited==t:
                    l.active=0
                    l=l.parent

        if verbose.find("2") <> -1:
            print "Candidates",cand
        

        # processing bt candidates
        for c in cand:
            if not c.l.active and not c.r.active:
                if c.interval[0].nearestt.lca(t)==t:
                    if verbose.find("2") <> -1:                    
                        print " cand ",c,"moved to (parentt)",t.parentt
                    t.parentt.bt.append(c)
                    #??? Nono (root T)?
                else:
                    #if k<0: print "Dupablada"    
                    c.interval[0].nearestt.bt.append(c)
                    if verbose.find("2") <> -1:                    
                        print " cand ",c,"moved to (>parentt)",c.interval[0].nearestt
                             
        if verbose.find("2") <> -1:        
            print "Remaining dupl:"   
        for d in dup:
            if d.active:
                if verbose.find("2") <> -1:                
                    print " ",d.num,d#,hasattr(d,'visited')
        if prevminscore>=0 and mescore>=prevminscore:
            if verbose.find("4") <> -1: 
                print "Skipping score ",mescore,"vs",prevminscore
                return -1
            return mescore 
    if verbose.find("3") <> -1: print "MEscore",mescore
    return mescore
            
        
     

ModLCA=1                # LCA model
ModGuigo=2              # GMS model Guigo et al.
ModPaszekGorecki=3      # PG model
ModFellows=4            # FHS model Fellows et al.

from itertools import combinations,chain


def ppscores(st):
    def speciestreescore(n):
        if n.leaf(): s=n.clusterleaf
        else: s="("+speciestreescore(n.r)+","+speciestreescore(n.l)+")"
        if n.scorek>0:
            s=s+" k=%d "%n.scorek
        return s
    total=sum(n.scorek for n in st.nodes)
    return speciestreescore(st.root)+" treename='ME=%d'"%total


def all_subsets(ss):
    return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))

def merfellows(gtrees,st):
    global verbose 
    lastopt=''
    if verbose.find("3") <> -1:
        print "@ Debug of Fellows model"
    specnodes=[ n for g in gtrees for n in g.root.nodes() if not n.leaf() and (n.lcamap!=n.l.lcamap and n.lcamap!=n.r.lcamap) ]
    if verbose.find("3") <> -1:
        # print "@ List of all speciation nodes before preprocessing",specnodes
        print "@ Number all speciation nodes before preprocessing",len(specnodes)  
    
    r=[]
    for s in specnodes:
        #if len(s.nodes())>3: 
        #print "SPEC:",s,s.nodes()
        hasdup=0
        for n in s.nodes():
            if not n.leaf() and (n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap):
                hasdup=1
                break
        if hasdup: r.append(s)
    specnodes=r



    subnodes=[]

    if verbose.find("3") <> -1: 
        print "Found ",len(specnodes),"speciation nodes"

    if verbose.find("5") <> -1: 
        print "@ List of all speciation nodes",specnodes

    maxdup=0
    alldupscore=0
    pompom=None
    pom=None
    pompomnum=0
    pomnum=0
    maxmaxdup=0
    change=False
    eq = False

    for i,g in enumerate(gtrees):
        maxh=maxc=-1
        change = False
        eq = False
        for n in g.root.nodes():
            if n.leaf():
                c=0
        
                while n:
                    if not n.leaf() and (n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap):
                        c=c+1                    
                    n=n.parent
            if c > maxdup: change = True
            if c == maxmaxdup: eq = True                    
            maxdup=max(maxdup,c)
            maxc=max(c,maxc)
        if change:
            maxmaxdup = maxdup
            pom = g
            pomnum = 0
        if eq and not change:
            pomnum = pomnum + 1            
        if alldupscore < g.height(): 
            pompom = g
            pompomnum = 0
        if alldupscore == g.height(): pompomnum = pompomnum + 1    
        alldupscore=max(alldupscore,g.height())
        if verbose.find("6") <> -1:      
            print i,maxh,maxc

    if verbose.find("1") <> -1: 
        print "Lower bound for ME score:",maxdup
        print "Lower bound tree:",pom
        print "Upper bound for ME score:",alldupscore
        print "Upper bound tree:",pompom
        if pom==pompom:
            if (pompomnum == 0) and (pomnum == 0):
                print "Unique tree G*. Unique tree which infers lower bound. Unique tree which infers upper bound."

    if verbose.find("8") <> -1: 
        print "Species tree"
        print str(st)
        for n in st.root.nodes():
            print "#%d %s\n"%(n.num,n)
        print "Speciation nodes after preprocessing from gene trees (number of node; number of lca-maped node in species tree; node)"
        for g in specnodes:
                print "#%d #%d %s\n"%(g.num,g.lcamap.num,g)

    if delnodes:
        print "@ List of all speciation nodes",specnodes
        print "@ Number of all speciation nodes",len(specnodes)
        print "Removing ",delnodes
        ll = delnodes.replace("[","").replace("]","").replace(",","") # strip('[],').split()
        for i,g in enumerate(gtrees):
            for n in g.nodes:
                if str(n.num) in ll.split(): 
                    subnodes.extend(n.nodes())
        print subnodes            
        autodups=[ n for n in specnodes if n in subnodes ]
        specnodes=[ n for n in specnodes if n not in subnodes ]
        print "@ List of all speciation nodes",specnodes
        print "@ Number of all speciation nodes",len(specnodes)


    if verbose.find("7") <> -1: 
        print "Found ",len(specnodes),"speciation nodes"
        if len(specnodes) > 20:
            quit()




    mescore=sum(len(gt.nodes) for gt in gtrees)+1
    first=1
    cnt=0
    for ss in all_subsets(specnodes):
    #for ss in [ specnodes ]:
        if verbose.find("3") <> -1:
            print "@ Subset of speciation nodes enabled for a change into duplication nodes",ss
        ss=list(ss)
        s = [n for n in specnodes if n not in ss ]
        #print len(s)  
        # zrob intervaly
        for gt in gtrees:
            for n in gt.nodes:
                n.interval=None
                if n.leaf(): continue
                if n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap:
                    n.interval=[n.lcamap, None]

        for n in s: n.interval=[n.lcamap,None]

        #HACK - given nodes deleted from SPEC, are automatically converted into duplications
        if delnodes:
            for n in autodups: n.interval=[n.lcamap,None]


        for gt in gtrees:
            expandintervals(gt,st)

        if verbose.find("6") <> -1:
            print "@ List of all intervals (spec->dup node intervals are extended)"        
            for gt in gtrees:            
                for n in gt.nodes:                
                    if n.interval:                    
                        print "Node ",n," interval ",n.interval     
        
        if first:
            m=mescore=mer(gtrees,st)
            first=0
            print "Current score",m
        else:
            m=meropt(gtrees,st,mescore)
        if m>-1: 
            if verbose.find("3") <> -1: print "%d. MEcurrent"%cnt,m,mescore
            if mescore>m:
                oldverbose,verbose=verbose,"1"
                print "*"*80
                print cnt, "REPEATED OPTIMAL COMP FOR ",m
                meropt(gtrees,st,-1)                
                lastopt=ppscores(st)
                print "*"*80                
                verbose=oldverbose
                print "Current score",m
            mescore=min(mescore,m)
        cnt=cnt+1
        if cnt%1000==0: print cnt,"variants processed; current min",mescore,"last score",m
        #if cnt==5000: quit()
    if printstreewithscores: print "&s",lastopt
    else: print "MEscore",mescore
    return mescore

    

def genLCAIntervals(gt,st):
    for n in gt.nodes:
        n.interval=None
        if n.leaf(): continue
        if n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap:
            n.interval=[n.lcamap,n.lcamap]  

def genPaszeGoreckiIntervals(gt,st):
    for n in gt.nodes:
        n.interval=None
        if n.leaf(): continue
        if n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap:
            n.interval=[n.lcamap, None]  

    expandintervals(gt,st)
    
def expandintervals(gt,st):
    for n in gt.nodes:
        if n.interval:        
            c=n.parent
            while c and c.interval: c=c.parent
            if not c: 
                n.interval[1]=st.root
                continue
            top=n.interval[0]
            while top.parent!=c.lcamap:
                top=top.parent
            n.interval[1]=top            


def genGMSIntervals(gt,st):
    for n in gt.nodes:
        n.interval=None
        if n.leaf(): continue
        if n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap:
            n.interval=[n.lcamap, None]  


    for n in gt.nodes:
        if not n.interval: continue
        if not n.parent:  # root
            n.interval[1]=st.root
            continue        

        if n.parent.interval and n.lcamap==n.parent.lcamap:
            n.interval[1]=n.interval[0]           
            continue
            
        top=n.lcamap
        

        while top.parent!=n.parent.lcamap:
            top=top.parent

        n.interval[1]=top            






def usage():
    print "Usage:"
    print " ********    MODELS OF GENERATING VALID MAPPINGS INTERVALS           ******** "
    print " -L             - lca-mapping model of valid mappings (least flexible)  "  
    print " -G             - GMS model of valid mappings"
    print " -P             - PG model of valid mappings"
    print " -F             - FHS model of valid mappings"
    print " ********    INPUT                                                   ********  "
    print "                - default argument is an input file (gene tree, species tree, intervals)"
    print " * note: please set model how to generate itervals and give both input species tree and gene tree *"
    print " -g TREE        - to define input gene tree"
    print " -s TREE        - to define input species treee"
    print " -p 'TREE TREE' - to give pair gene tree, species tree"
    print " -m FILE        - defines a set of gene trees"
    print " -j NUMTREE,LEN - skip first NUMTREE trees and take next LEN trees; missing LEN=all trees"    
    print " -t             - for each gene tree print its height"
    print " -k             - print ME score and a species tree with partial ME scores at corresponding nodes"
    print " -i             - set intervals from file"    
    print " ********    OUTPUT                                                  ******* "
    print " -o             - write to default output file"
    print " -O FILE        - write to user defined output file"
    print " -v STRING      - verbose mode; (default) print ME score only; debug 0 - input/output; 1 - upper, lower bounds for ME score; 2 - detailed Alg.2 debug; 3 - Alg.3 debug; 4 - skipped scores;5 - print nodes from SPEC; 6 - more Alg.3 details; 7 - print size of SPEC and stops if > 22; 8 - print node number for species tree and speciations from SPEC; where SPEC is a set of all lca-speciation nodes (from input gene trees) which are ancestor of some duplication."
    print " -d STRING      - list of node numbers to delete from SPEC and convert into duplications"
    

def readgtreefile(gtreefile):

    gtrees=[]
    for g in open(gtreefile).readlines():
         if g and g[0]!="#":
            gtrees.append(Tree(str2tree(g.strip())))

    return gtrees




def main():
     # parse command line options

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ktm:ep:j:g:s:FPGLioO:hv:d:", ["help", "output="])
    except getopt.GetoptError as err:
        print str(err) 
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        print "Usage: [-e] [-p] [-r] inputfile"
        sys.exit(1)

    model=3 # default PG model
    stree = None
    gtrees = []
    gt=st=None
    global verbose,printstreewithscores,delnodes

    outputfile=None
    runmer=1 # default run
    gtreefile=None
    firstgtrees=-1
    numtrees=-1
    ppheight=0
    delnodes=""
    printstreewithscores=0
    for o, a in opts:
        if o == "-s":
            st = Tree(str2tree(a)) 
        elif o == "-v":
            verbose=str(a)
        elif o == "-d":
            delnodes=str(a)
        elif o == "-p":
            gtree,stree=a.strip().split(" ")
            gtrees.append(Tree(str2tree(gtree)))
            st=Tree(str2tree(stree))
            
        elif o == "-g":
            gtrees.append(Tree(str2tree(a)))
            
        elif o == "-m":
            gtreefile=a 
        elif o == "-i":
            model=0


        elif o == "-t":
            ppheight=1
        elif o == "-G":
            model=ModGuigo
        elif o == "-P":
            model=ModPaszekGorecki

        elif o == "-R":
            model=ModRandom
        elif o == "-L":
            model=ModLCA
        elif o == "-F": 
            model=ModFellows
        elif o =="-o":
            outputfile="outmer.txt"
        elif o =="-O":
            outputfile=a

        elif o=="-j":
            if "," in a: firstgtrees,numtrees=int(a.split(",")[0]),int(a.split(",")[1])
            else: firstgtrees,numtrees=int(a),-1


        elif o=="-k":
            printstreewithscores=1
            
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"

    if gtreefile: 
        gtrees.extend(readgtreefile(gtreefile))
    
    for l in args:
        gtrees,st=readintervalfile(l)

    if ppheight:
        for gt in gtrees:
            print gt.height(),gt
        if not st: sys.exit(0)

    if not st or not gtrees:
        print "Both trees have to be defined"
        usage()
        sys.exit(3)

    if firstgtrees>=0:
        if numtrees>0:
            gtrees=gtrees[firstgtrees:firstgtrees+numtrees]
        else: 
            gtrees=gtrees[firstgtrees:]




    if model:
        for gt in gtrees:   
            gt.setlcamapping(st)

    for gt in gtrees:  
        if model==ModPaszekGorecki:
            genPaszeGoreckiIntervals(gt,st)
        elif model==ModGuigo:
            genGMSIntervals(gt,st)
        elif model==ModLCA:
            genLCAIntervals(gt,st)


    if model==ModFellows:
        merfellows(gtrees,st)
        sys.exit(0)

    if outputfile: 
        savegsi(gtrees,st,outputfile)

    if verbose == 3:     
        for i,gt in enumerate(gtrees):        
            for g in gt.root.nodes():            
                if g.interval:                
                    print "Interval for %d in tree%d:"%(g.num,i),g,g.interval

        
    if runmer:
        me=mer(gtrees,st)
        if printstreewithscores: print "&s",ppscores(st)
        else: print "MEscore",me


    return 0
        

if __name__ == "__main__":
    main()
