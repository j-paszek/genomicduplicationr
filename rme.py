#!/usr/bin/python

import sys
from treeop import Tree, str2tree,  node2label
import itertools


verbose = ""
printstreewithscores = 0

def readgs(filename):

    t=[ l.strip() for l in open(filename,"r").readlines() if len(l.strip()) and l.strip()[0]!='#' ]
    if verbose == 3:
        print(t)
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
    print("File %s saved"%outputfile)


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

    if verbose == 1:    
        print("Duplication Nodes",dup)

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


    if verbose == 1:
        print("Nodes of T")
        for t in tdup: 
            print(t)
            print("   Bt=",t.bt,[ n.num for n in t.bt])
            print("   Gammat=",t.duproots)

    def travduproots(r):
        if r.leaf(): return 1
        
        h1=h2=0
        if r.l.active: h1=travduproots(r.l)
        if r.r.active: h2=travduproots(r.r)
        return 1+max(h1,h2)

    if verbose == 1:    
        print("Dup Leaves",ldup)

    mescore=0


    # Main loop
    for t in tdup:
        if verbose == 1:        
            print("="*80)
            print("Processing",t)
            print("   Bt=",t.bt)
            print("   Gammat=",t.duproots)

        
        k=-1        
        for r in t.duproots:
            if r.active:
                k=max(k,travduproots(r))
        if verbose == 1:        
            print("   k=",k)

        if k>0: 
            mescore+=k
            t.scorek=k
        else: t.scorek=0

    
        cand=[]
        if k<0:
            if verbose == 1:            
                for l in t.bt:
                    print("Moving from orhpan dup leaf",l,l.active)
            
            t.parentt.bt.extend(t.bt)

        else:
            #processing leaves from Bt, k>0 case 
            
            for l in t.bt:
                if verbose == 1:                 
                    print("Removing from dup leaf",l.num,l,l.active)
                curk=1
                newcand=None
                while curk<=k and not newcand:
                    l.visited=t
                    if verbose == 1:                    
                        print("  Marking",l.num,l)
                    if not l.parent: break # gene tree root 

                    l.h=curk
                    sib=l.sibling()
                    l=l.parent
                    if verbose == 1:
                        print("  checking",l,l.num if l else "") # l.interval
                    if not l.active:
                        if verbose == 1: 
                            print("  not active")
                        break
                    if verbose == 1:
                        print(l.interval)
                    
                    if l.interval[0].nearestt.lca(t)!=t or curk==k:
                        newcand=l
                        if verbose == 1: print(" cand. found",l)
                        # stop
                    elif not sib.active:
                        curk+=1
                        if verbose == 1: print("  sib is not active - cont",sib)
                    elif sib.visited==t:
                        curk=max(curk,sib.h)+1                
                    else:
                        if verbose == 1:                        
                            print("  sib not visited yet",sib)
                        break

                if newcand: cand.append(newcand)

            # Clean active nodes
            for l in t.bt:
                while l and l.visited==t:
                    l.active=0
                    l=l.parent

        if verbose == 1:
            print("Candidates",cand)
        

        # processing bt candidates
        for c in cand:
            if not c.l.active and not c.r.active:
                if c.interval[0].nearestt.lca(t)==t:
                    if verbose == 1:                    
                        print(" cand ",c,"moved to (parentt)",t.parentt)
                    t.parentt.bt.append(c)
                    #??? Nono (root T)?
                else:
                    #if k<0: print "Dupablada"    
                    c.interval[0].nearestt.bt.append(c)
                    if verbose == 1:                    
                        print(" cand ",c,"moved to (>parentt)",c.interval[0].nearestt)
                             
        if verbose == 1:        
            print("Remaining dupl:")
        for d in dup:
            if d.active:
                if verbose == 1:                
                    print(" ",d.num,d)#,hasattr(d,'visited')
        if prevminscore>=0 and mescore>=prevminscore:
            if verbose == 4: 
                print("Skipping score ",mescore,"vs",prevminscore)
                return -1
            return mescore 
    if verbose == 2: print("MEscore",mescore)
    return mescore
            
        
     


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
    if verbose ==4:
        print("@ Debug of Fellows model")
    specnodes=[ n for g in gtrees for n in g.root.nodes() if not n.leaf() and (n.lcamap!=n.l.lcamap and n.lcamap!=n.r.lcamap) ]
    if verbose == 4:
        # print "@ List of all speciation nodes before preprocessing",specnodes
        print("@ Number all speciation nodes before preprocessing",len(specnodes))
    
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

    if verbose ==4:
        print("@ List of all speciation nodes",specnodes)
        print("@ Number of all speciation nodes",len(specnodes))
        #quit()



    subnodes=[]
    autodups=[]
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

    if verbose ==2: 
        print("Found ",len(specnodes),"speciation nodes")
        print("@ List of all speciation nodes",specnodes)

    maxdup=0
    alldupscore=0
    pompom=None
    pom=None

    for i,g in enumerate(gtrees):
        maxh=maxc=-1
        for n in g.root.nodes():
            if n.leaf():
                c=0
        
                while n:
                    if not n.leaf() and (n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap):
                        c=c+1                    
                    n=n.parent
            if c > maxdup: pom = g
            maxdup=max(maxdup,c)
            maxc=max(c,maxc)
        if alldupscore < g.height(): pompom = g
        alldupscore=max(alldupscore,g.height())
        if verbose ==6:      
            print(i,maxh,maxc)

    if verbose ==2: 
        print("Bottom limit for ME score:",maxdup)
        print("Bottom limit tree:",pom)
        print("AllDupScore:",alldupscore)
        print("AllDup tree:",pompom)
        #if pom==pompom:
        #    print " Its the same tree !!!"
        # quit()

    if verbose ==6:
        #mer.py -F -j 1257,1 -v6 -e in.txt
        #w treefam

        if len(gtrees)>10: 
            print("too many trees")
            sys.exit(-1)  # stop    

        def ptree(t):
            if t.leaf(): 
                t.lf=1

                return
            if t.lcamap==t.l.lcamap or t.lcamap==t.r.lcamap:
                t.lf=0
                t.lcamap.ds|=1                
            elif t in specnodes:
                t.lcamap.ds|=2 
                t.lf=0 
            else: t.lf=1
            if not t.lf:
                ptree(t.l)
                ptree(t.r)

        def ppniceold(t,isst): # up to two marks in S
            
            l=''
            dp=0
            if not isst:

                if t.leaf(): 
                    return "lw=0.2 s=\"%s\" :0.0"%t.clusterleaf
                if (t.lcamap==t.l.lcamap or t.lcamap==t.r.lcamap):
                    l="mark=d%d nn=\"%d\""%(t.lcamap.num,t.num)
                    dp=1
                if t in specnodes:
                    l="mark=s%d nn=\"%d\""%(t.lcamap.num,t.num)
                                
                if not l:
                    #return "mark=94 :0.3"#%(t.num,t.num)
                    return ":0.0"#+("mark=100" if t.lf else "")
                    l="mark=s%d "%t.lcamap.num
                    t.lcamap.ds|=2        
                    return l+"s=''" #\"%d\""%t.lcamap.num
            
                #x=" s=\"!R!%d \""%t.lcamap.num
                x=''
            
                if t.lcamap.leaf() and dp==1:
                    if t.r.stheight>t.l.stheight:
                        return "("+ppniceold(t.r,isst)+")"+l
                    return "("+ppniceold(t.l,isst)+")"+l
            else:                 
                if t.ds==3: x=" mark=[d%d,s%d] "%(t.num,t.num)
                elif t.ds==2: x=" mark=[s%d] "%t.num
                elif t.ds==1: x=" mark=[d%d] "%t.num
                else: x=''
                x=x+" nn=\"%d\""%t.num
                
                if t.leaf(): return x+" s=\"%s\""%t.clusterleaf
                #else: x=x+" s=\"!R!%d \""%(t.num)
            if t.l.lf:
                return "("+ppniceold(t.r,isst)+") "+l+x
            elif t.r.lf: 
                return "("+ppniceold(t.l,isst)+") "+l+x
            return "("+ppniceold(t.l,isst)+","+ppniceold(t.r,isst)+") "+l+x



        def ppnice(t,isst): # single mark in S
            
            l=''
            dp=0
            maxmark=29*3
            if not isst:

                if t.leaf(): 
                    return "lw=0.2 s=\"%s\" :0.0"%t.clusterleaf
                if (t.lcamap==t.l.lcamap or t.lcamap==t.r.lcamap):
                    l="mark=d%d nn=\"%d\""%(t.lcamap.num,t.num)
                    t.lcamap.ds|=1
                    dp=1
                if t in specnodes:
                    l="mark=s%d nn=\"%d\""%(t.lcamap.num,t.num)
                    t.lcamap.ds|=2
            
                if not l:
                    #return "mark=94 :0.3"#%(t.num,t.num)
                    return ":0.0"#+("mark=100" if t.lf else "")
                    
                x=''
            
                if t.lcamap.leaf() and dp==1:
                    if t.r.stheight>t.l.stheight:
                        return "("+ppnice(t.r,isst)+")"+l
                    return "("+ppnice(t.l,isst)+")"+l
            else:                 
                if t.ds==3: x=" mark=[d%d,s%d] "%(t.num,t.num)
                elif t.ds==2: x=" mark=s%d "%(t.num)
                elif t.ds==1: x=" mark=d%d "%t.num
                else: x=''
                x=x+" nn=\"%d\""%t.num
                
                if t.leaf(): return x+" s=\"%s\""%t.clusterleaf
                #else: x=x+" s=\"!R!%d \""%(t.num)
            if t.l.lf:
                return "("+ppnice(t.r,isst)+") "+l+x
            elif t.r.lf: 
                return "("+ppnice(t.l,isst)+") "+l+x
            return "("+ppnice(t.l,isst)+","+ppnice(t.r,isst)+") "+l+x
        f=open("dup.gse","w")
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
            n.ds=0
            n.lf=0
            n.sp=n.dp=0
        for g in gtrees:
            for n in g.nodes: n.lf=0       
            ptree(g.root)
            cnt=0
            for n in st.nodes:
                if n.ds or n.sp: 
                    n.marknum=cnt
                    cnt+=1
            f.write("&g %s\n"%ppnice(g.root,0))            
        f.write("&s %s\n"%ppnice(st.root,1))
        cnt=0
        maxmark=29*3
        for n in st.nodes: 
            if n.ds&1: 
                f.write("d%d=%d\n"%(n.num,n.marknum))
            if n.ds&2: 
                f.write("s%d=%d\n"%(n.num,n.marknum+maxmark))
                

        # print "Dupnodesmap:",cnt,"Specnodesmap",scnt-23*3

        f.close()


        print("File dup.gse saved")

        quit() 
    mescore=sum(len(gt.nodes) for gt in gtrees)+1
    first=1
    cnt=0
    for ss in all_subsets(specnodes):
    #for ss in [ specnodes ]:
        if verbose ==4:
            print("@ Subset of speciation nodes enabled for a change into duplication nodes",ss)
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

        #HACK
        for n in autodups: n.interval=[n.lcamap,None]


        for gt in gtrees:
            expandintervals(gt,st)

        if verbose == 4 :
            print("@ List of all intervals (spec->dup node intervals are extended)")
            for gt in gtrees:            
                for n in gt.nodes:                
                    if n.interval:                    
                        print("Node ",n," interval ",n.interval)
        
        if first:
            m=mescore=mer(gtrees,st)
            first=0
            print("Current score",m)
        else:
            m=meropt(gtrees,st,mescore)
        if m>-1: 
            if verbose ==4: print("%d. MEcurrent"%cnt,m,mescore)
            if mescore>m:
                oldverbose,verbose=verbose,1
                print("*"*80)
                print(cnt, "REPEATED OPTIMAL COMP FOR ",m)
                meropt(gtrees,st,-1)                
                lastopt=ppscores(st)
                print("*"*80)
                verbose=oldverbose
                print("Current score",m)
            mescore=min(mescore,m)
        cnt=cnt+1
        if cnt%1000==0: print(cnt,"variants processed; current min",mescore,"last score",m)
        #if cnt==5000: quit()
    if printstreewithscores: print("&s",lastopt)
    else: print("MEscore",mescore)
    return mescore

    

def genLCAIntervals(gt,st):
    for n in gt.nodes:
        n.interval=None
        if n.leaf(): continue
        if n.lcamap==n.l.lcamap or n.lcamap==n.r.lcamap:
            n.interval=[n.lcamap,n.lcamap]  

def genPaszekGoreckiIntervals(gt,st):
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
        #if n.num==19: print "ASDADA",n
        if not n.parent:  # root
            n.interval[1]=st.root
            continue        

        if n.parent.interval and n.lcamap==n.parent.lcamap:
            n.interval[1]=n.interval[0]           
            continue
            
        top=n.lcamap
        

        while top.parent!=n.parent.lcamap:
            #if n.num==19: 
            #    print "-----",top,n.lcamap,n.parent.lcamap
            top=top.parent
        #if n.num==19: 
            #print "--LAST",top,n.lcamap,n.parent.lcamap
        n.interval[1]=top            

def ppgsevoloutput(gtrees,st):
    f=open("mer.gse","w")
    f.write("#!/home/gorecki/Dropbox/projects/htree/gsevol.py -x\n")
    for i,gt in enumerate(gtrees):
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


"""%(i,gt,st))

    f.close()

    
def readgtreefile(gtreefile):

    gtrees=[]
    for g in open(gtreefile).readlines():
         if g and g[0]!="#":
            gtrees.append(Tree(str2tree(g.strip())))

    return gtrees


