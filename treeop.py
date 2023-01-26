import sys

def leaf(s): return type(s)==str
def comparable(a,b): set(a).intersection(set(b))
def num(c): return ordext(c)


def height(t):
    if leaf(t): return 0
    return 1+max(height(t[0]),height(t[1]))
           
def theight(t,l):

    def _th(t,l,h):
        if leaf(t):
            if l==t: return h
            return 0

        mx=_th(t[0],l,h+1)
        if mx>0: return mx
        return _th(t[1],l,h+1)

    return _th(t,l,0)




lfalf=[ chr(i) for i in xrange(ord('a'),ord('z')) ]+[ chr(i) for i in xrange(ord('A'),ord('Z')) ]
#lfalf="ABCDEF" #abcdyzghjiklmopqr123456"

lfalf=[ chr(i) for i in xrange(ord('a'),ord('z')) ]+[ chr(i) for i in xrange(ord('A'),ord('Z')) ]
def ordext(c): 
    return lfalf.index(c)


# krotka -> napis
def pt(s): return str(s).replace("'",'').replace(" ",'')

# napis -> krotka
def str2tree(s):

    def _st(s,p):
        if not s:
            raise Exception("String too short: <%s>"%s)
        if s[p]=='(':
            t,p=_st(s,p+1)
            if s[p]!=',': raise Exception("comma expected")
            t2,p=_st(s,p+1)
            if s[p]!=')': raise Exception(") expected")
            return ((t,t2),p+1)
        r=''
        while s[p].isalnum():
            r+=s[p]
            p=p+1
        return (r,p)

    if not s.strip():
        print "Warning: empty string"
        return []
    return _st(s,0)[0]


def node2label(n):
    return str(n).replace("("," ").replace(")"," ").replace(","," ")


lfalf=[ chr(i) for i in xrange(ord('a'),ord('z')) ]+[ chr(i) for i in xrange(ord('A'),ord('Z')) ]


def getclusters(s):
    def _gc(s,d):
        d[s]=cluster(s)
        if not leaf(s):
            _gc(s[0],d)
            _gc(s[1],d)
    d={}
    _gc(s,d)
    return d

class Node:
    def __init__(self, tup, par):
        self.src=tup            
        self.parent=par

        if self.parent:
            self.height=self.parent.height+1
        else: self.height=0
        self.stheight=height(tup)
        if self.leaf():
            self.cluster=frozenset([tup])
            self.clusterleaf=tup
        else:
            self.l=Node(tup[0],self)
            self.r=Node(tup[1],self)
            self.cluster=self.r.cluster | self.l.cluster
            self.clusterleaf=None

        self.interval=None

    def leaf(self):
        return type(self.src)==str

    def setinterval(self,s,distornode):
        if type(distornode)==int:
            self.interval=[s,s.ancestor(distornode)]
        else: self.interval=[s,distornode]

    def ancestor(self,dist):
        if dist==0 or not self.parent: return self
        return self.parent.ancestor(dist-1)

    def __str__(self):
        if self.leaf(): return self.clusterleaf
        return "("+self.l.__str__()+","+self.r.__str__()+")"

    def __repr__(self):
        return self.__str__()


    def nodes(self): # postorder
        if self.leaf(): return [ self ]
        return  self.l.nodes()+ self.r.nodes() + [self]

    def leaves(self):
        if self.leaf(): return [ self ]
        return  self.l.leaves()+ self.r.leaves()
    
    def sibling(self):
        if not self.parent: return None
        if self.parent.l==self: return self.parent.r
        return self.parent.l
    
   # x = self or x below sel.
    def geq(self,x):
        while x:
            if x==self: return True
            x=x.parent
        return False

    def leq(self,x):
        return x.geq(self)

    def comparable(self,x):
        return self.geq(x) or self.leq(x)

    def lca(self,y):
        a,b=self,y
        if a.height>b.height: a,b=b,a
        # a.h <= b.h        
        while b.height!=a.height: b=b.parent

        # a.h == b.h        
        while True:
            if a==b: return a
            a=a.parent
            b=b.parent

    def findnode(self,cluster):

        #print "FD",self,cluster,self.src

        if self.cluster==cluster: return self
        if self.leaf(): return None
        r=self.l.findnode(cluster)
        if r: return r
        return self.r.findnode(cluster)



class Tree:
    def __init__(self,tup):
        self.srclist=None
        #print type(tup),tup
        if type(tup)==list:
            self.srclist=tup
            tup=self.srclist[0]
            for t in self.srclist[1:]: tup=(tup,t)

        self.root=Node(tup,None)
        self.nodes=self.root.nodes()
        for i,n in enumerate(self.nodes): 
            n.num=i
            n.artificial=False
        self.src=tup

        if self.srclist:
            l=self.root
            for i in xrange(len(self.srclist)-1):
                #print "!",l
                l.artificial=True
                l=l.l
        
    def leaves(self):
        return self.root.leaves()


    def findnodeplus(self,cluster):

        if len(cluster)>2 and cluster[-2]=="+":
            up=int(cluster[-1])
            cluster=cluster[:-2]
        elif len(cluster)>1 and cluster[-1][0]=="+":
            up=int(cluster[-1])
            cluster=cluster[:-1]
        else:
            up=0

        s=self.root.findnode(frozenset(cluster))
        if not s: return s

        while up and s.parent:
            s=s.parent
            up=up-1
        return s

    def __str__(self):
        return self.root.__str__()

  
    def setlcamapping(self,st):
        for l in self.leaves():
            clu=l.clusterleaf
            l.lcamap=None
            while clu and not l.lcamap:
                l.lcamap=st.root.findnode(frozenset([clu]))
                clu=clu[:-1]
            if not l.lcamap:
                raise Exception("Lca mapping not found for ",l)

        for l in self.nodes:
            if not l.leaf():
                l.lcamap=l.r.lcamap.lca(l.l.lcamap)
            #print "LCA",l,l.lcamap

        
    def __repr__(self):
        return self.root.__repr__()

    
    def lcacluster(self,cluster):
        c=set(cluster)
        for n in self.nodes:
            if c.issubset(set(n.cluster)): return n
     
    
    def clusters(self):
        return set(n.cluster for n in self.nodes)


    def height(self):
        return self.root.stheight



        
