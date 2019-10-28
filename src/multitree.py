
import random
import itertools
import sys

lfalf=[ chr(i) for i in range(ord('a'),ord('z')) ]+[ chr(i) for i in range(ord('A'),ord('Z')) ]


gencomment="""
Maddisson formula:

DC(G,S) is the sum of k(v) over all nodes of G except the root,
where
   k(v)=||M(v),M(parent(v))||-1,
   M:G->S is the lca-mapping
and
   ||v,w|| is the number of edges on the path connecting v and w in S.

The costs from the MAX DC paper should be adjusted by adding 2*|G|-2.
"""


def cluster(s): 
    if type(s)==str: return set([s])
    return cluster(s[0]).union(cluster(s[1]))


# tuple -> string 
def pt(s): return str(s).replace("'",'').replace(" ",'')

# string -> tuple 
def str2tree(s):
    def _st(s):
       
        s=s.strip()
        if not s:
            raise Exception("String too short: <%s>"%s)
        if s[0]=='(':
            t1,s=_st(s[1:])
            lst=(t1,)
            while s[0]==',':                
                t1,s=_st(s[1:])
                lst=lst+(t1,)
            if s[0]!=')': raise Exception(") expected")
            return (lst,s[1:])
        
        lab=''
        while s and s[0].isalnum():
            lab=lab+s[0]
            s=s[1:]
        if not lab:
            print("Label expected in tree string")
            sys.exit(1)
        return (lab,s)

    if not s.strip():
        print("Warning: empty string")
        return []

    return _st(s)[0]


def randtree(n,chnum=0):
    l=lfalf[:]
    c=0
    while n>len(l):
        l.extend(lab+"%d"%c for lab in lfalf)
        c+=1
    return _randtree(l[:n],chnum)
           

def randmultitree(t,mprob):
    if type(t)==str: t=str2tree(t)
    return _randmultitree(t,mprob)

def _randmultitree(t,mprob):
    if type(t)==str: return t    
    t=tuple([_randmultitree(c,mprob) for c in t])
    res=[]
    for c in t:
        if c!=str and random.random()<mprob: 
            # flattenif             
            res.extend(list(c)) 
            continue
        res.append(c)            
    return tuple(res)

def _randtree(lst,chnum=0):
    if len(lst)==1: return lst[0]
    if chnum<=0:
        dp=random.randint(1,len(lst)-1)
        return (_randtree(lst[:dp]),_randtree(lst[dp:]))

    lf=lst[:]
    ch=[]

    if chnum:
        cher=[]
        for i in range(chnum):
            cher.append((lf.pop(random.randint(0,len(lf)-1)),lf.pop(random.randint(0,len(lf)-1))))
        lf.extend(cher)        

    while len(lf)>1:        
        curl1=lf.pop(random.randint(0,len(lf)-1))
        curl2=lf.pop(random.randint(0,len(lf)-1))            
        if chnum and type(curl1)==str==type(curl2): 
            lf.append(curl1)
            lf.append(curl2)
            continue 
        lf.append((curl1,curl2))

    return lf[0]


class Tree:
    def __init__(self,tup):
        if type(tup)==str: tup=str2tree(tup)
        self.root=Node(tup,None)
        self.nodes=self.root.nodes()
        for i,n in enumerate(self.nodes): 
            n.num=len(self.nodes)-i-1
            n.nodeid=i
        self.src=tup
        self.n=len(self.root.cluster)

        # print  tup,self.n
        # for n in self.nodes:
        #     print n,n.cluster
        # print "LEAV",self.leaves()

    def weight(self):
        return sum( n.depth for n in self.leaves() )


    def leaves(self):
        return self.root.leaves()

    def __str__(self):
        return self.root.__str__()


    def lcacluster(self,cluster):
        c=set(cluster)
        for n in self.nodes:
            if c.issubset(set(n.cluster)): return n


    def inferinternalmaps(self):
        for g in self.nodes:
            if g.c:
                newmap = g.c[0].map
                for c in g.c:
                    newmap = newmap.lca(c.map)
                g.map = newmap
                #g.map=reduce(lambda s,g: s.lca(g.map),g.c,g.c[0].map)

    def setlcamapping(self,st):

        #def setlcamapping(gt,st,**kwargs):
        stleaves=st.leaves()
        #print "STL",stleaves
        for n in self.leaves():
            #print "LCA",n,n.cluster
            stl=[s for s in stleaves if s.cluster==n.cluster ]
            if len(stl)!=1:
                raise Exception("Ambiguous mapping: %d candidates. Label: %s" % (len(stl),n.cluster ))
            n.map=stl[0]            

        self.inferinternalmaps()        

    def cut(self,lf):
        return self.root.cut(lf)
    def cutedges(self):
        return self.root.cutedges()
            

    # Naive alg. 
    def genrank(self,s):


        for n in s.nodes: n.rank=n.height+1
        self.setlcamapping(s)                    
        for g in self.nodes: 

            if not g.c: continue
            if len(g.c)!=2: 
                raise Exception("Binary gene tree expected. Found %s"%g)

            r=100000
            for i in g.c[0].leaves():
                for j in g.c[1].leaves():
                    r=min(r,i.map.lca(j.map).rank)
            g.esr=r
            #print g,r
            #_setlabel(g,'esr',"%s"%r)
            #_setlabel(g,'lcarank',"%s"%g.lcamapping.rank)

    # O(d|G|logd +|S|)
    def genrankbyalg(self,s,verboselevel=0):
        sl=s.leaves()
        gl=self.leaves()
        for c in sl: c.glist=[]
        for i,c in enumerate(gl): 
            c.id=i+1
            c.map.glist.append(c)
        
        for n in self.nodes: n.esrbyalg=0
        snodes=sorted(s.nodes[:],lambda x,y: x.rank-y.rank) # can be directly traversed

        def np(Lambda):
            return " ".join( "["+",".join("%s%d"%(k,k.id) for k in l)+"]" for l in Lambda)

        def merge(n,lst):
            res=[]            
            
            if verboselevel & 8: 
                print("==== MERGE %s rank=%d Lambda=%s"%(n,n.rank,np(lst)))
            while lst:                
                minleafid=min(lst[i][0].id for i in range(len(lst)))
                for i in range(len(lst)):
                    if lst[i][0].id==minleafid:
                        lst[i].pop(0)
                        if not lst[i]: lst.pop(i)
                        break
                minleaf=gl[minleafid-1]
                if res:                           
                    lcagnode=res[-1].lca(minleaf)             

                    if verboselevel & 8: 
                        print( " Checking %s%d %s%d g=%s g.rank=%d:"%(res[-1],res[-1].id,
                            minleaf,minleaf.id,lcagnode,lcagnode.esrbyalg),)
                    if not lcagnode.esrbyalg: 
                        if verboselevel & 8: 
                            print("  Rank SET!",n.rank)

                        lcagnode.esrbyalg=n.rank
                        lcagnode.minr=[res[-1],minleaf]                                                  

                    else: 
                        if verboselevel & 8: print(" Rank Ignored :(")
                res.append(minleaf)

            if verboselevel & 8: 
                print("MergeResult",np([res]))
            return res

        
        for n in snodes:
            if n.leaf(): n.glist=merge(n,[n.glist])
            else: n.glist=merge(n,[ c.glist for c in n.c])

        # Check correctness with naive alg
        for n in self.nodes:
            if not n.leaf():
                #print n.esr,n.esrbyalg,n
                if n.esr!=n.esrbyalg:
                    print("ERR")
                    sys.exit(-1)
    


    # O(d|G|+|S|) algorithm 
    def genrankbyalg2(self,s,verboselevel=0):
        sl=s.leaves()
        glprefix=self.leaves()
        for c in sl: c.glist=[]
        for i,c in enumerate(glprefix): 
            c.id=i+1
            c.smap=c.map # set to its label mapping
        
        for n in self.nodes: n.esrbyalg2=0 # init  
        for n in s.nodes: n.lastgleaf=0 # init
        d=max(n.rank for n in s.nodes)
        for r in range(1,d+1):
            if verboselevel & 8: print("="*30,"rank=%d"%r)
            for v in glprefix:
                if verboselevel & 8: print(v,v.smap,v.smap.rank)
                if v.smap.rank!=r: continue
                if v.smap.lastgleaf: 
                    l=v.smap.lastgleaf
                    lcagnode=v.smap.lastgleaf.lca(v)  
                    
                    if verboselevel & 8: 
                        print("Checking %s%d %s%d g=%s g.rank=%d:"%(l,l.id,
                            v,v.id,lcagnode,lcagnode.esrbyalg2),)
                    
                    if not lcagnode.esrbyalg2:
                        if verboselevel & 8: 
                            print("  Rank SET!",r)
                        lcagnode.esrbyalg2=r # set rank
                        lcagnode.minr=[v.smap.lastgleaf,v]
                v.smap.lastgleaf=v
                v.smap=v.smap.parent # climb

        # Check correctness with naive alg
        for n in self.nodes:
            if not n.leaf():
                if n.esr!=n.esrbyalg2:
                    print(n.esr,n.esrbyalg2,n)
                    print("ERR2")
                    sys.exit(-1)      

    def ppgse(self):
        return self.root.ppgse()

    def dupcost(self,st):
        self.setlcamapping(st)
        c=0
        for n in self.nodes:            
            if n.leaf(): continue
            if n.c[0].map==n.map or n.c[1].map==n.map: c+=1
        return c
    
    def losscost(self,st):
        #tylko dla binarnych
        return self.dccost(st)+2*self.dupcost(st)

    def losscost2(self,st):
        # root has no losses
        # single loss in multifurcations        
        # duplications are not allowed (not checked here)
        self.setlcamapping(st)
        c=0
        for n in self.nodes:
            if n==self.root: continue        
            c+=n.map.depth-n.parent.map.depth-1 # single losses inside the path
            # check bottom node
            if not n.leaf() and len(n.map.c)>2: c+=1
        return c






    def dccost(self,st):
        self.setlcamapping(st)        
        c=0
        for n in self.nodes:                
            if n!=self.root: 
                c+=n.map.depth-n.parent.map.depth-1
        return c
            

class Node:
    def __init__(self, tup, par):
        self.src=tup
        self.cluster=cluster(tup)
        self.parent=par
        if self.parent:
            self.depth=self.parent.depth+1
        else: 
            self.depth=0

        if type(tup)==str:
            self.height=0
            self.c=None
            # leaf
        else:
            self.c = [ Node(t,self) for t in tup ]            
            self.height=max(c.height for c in self.c)+1



    def siblinggen(self):
        if self.parent:
            for i in self.parent.c: 
                if i!=self: yield i

    def siblingclustern(self):
        c=set([])
        for i in self.siblinggen(): c.union(i.cluster())
        return c
            

    def _cutedges(self):        
        if self.leaf(): 
            return self.label(),[]
        t1,f1=self.c[0]._cutedges()
        t2,f2=self.c[1]._cutedges()
        if self._cut==0:
            return (t1,t2),f1+f2
        if self._cut==1:
            return t2,f1+f2+[t1]
        return t1,[t2]+f1+f2

    def cutedges(self):
        t,f=self._cutedges()
        return [t]+f


    def cut(self,lf):
        if self.leaf():
            if self.label() in lf:
                return self.label()
            return None
        lst=[ c.cut(lf) for c in self.c ]
        lst=[ c for c in lst if c ]
        if len(lst)>1: 
            return "("+",".join(lst)+")"
        if len(lst)==1: return lst[0]
        return None
            

    def leaf(self):
        return not self.c       

    def label(self):
        return list(self.cluster)[0]

    def smp(self):
        return "".join(c.label() for c in self.leaves())

    def smpd(self):
        if self.leaf(): return self.label()        
        return "|".join(c.smp() for c in self.c)
        

    def __str__(self):
        if self.leaf(): return list(self.cluster)[0]
        return "("+",".join( c.__str__() for c in self.c)+")"

    def __repr__(self):
        return self.__str__()

    def nodes(self): # postfix
        if self.leaf(): return [ self ]
        # children first
        return  sum((c.nodes() for c in self.c),[])+[self]

    def leaves(self):
        if self.leaf(): return [ self ]
        return  sum((c.leaves() for c in self.c),[])
    

    def comparable(self,x):
        return self.geq(x) or self.leq(x)

    def leq(self,x):
        return x.geq(self)

    # x = self or x below sel.
    def geq(self,x):
        while x:
            if x==self: return True
            x=x.parent
        return False

    def lca(self,y):
        a,b=self,y
        if a.depth>b.depth: a,b=b,a
        # a.h <= b.h        
        while b.depth!=a.depth: b=b.parent

        # a.h == b.h        
        while True:
            if a==b: return a
            a=a.parent
            b=b.parent

    def ppgse(self): 
        s=''
        if hasattr(self,'esr'): s+=" esr=%d"%self.esr
        if hasattr(self,'rank') and not self.leaf(): s+=" rank=%s"%self.rank        
        #if hasattr(self,'minr'): s+=" minr='%s_%d-%s_%d' "%(self.minr[0],self.minr[0].id,self.minr[1],self.minr[1].id)
        if hasattr(self,'minr'): s+=" minr='%s_{%d}%s_{%d}'"%(self.minr[0],self.minr[0].id,self.minr[1],self.minr[1].id)
        if self.leaf(): 
            if hasattr(self,'id'): return self.label()+s+" leaflabel='%s_{%d}'"%(self.label(),self.id)
            return self.label()+s
        return "("+",".join(c.ppgse() for c in self.c)+")"+s


