'''
MPS test.
'''

__all__=['test_mps']

from numpy import *
from HamiltonianPy.Math.Tensor import *
from HamiltonianPy.DMRG.MPS import *
from copy import copy,deepcopy

def test_mps():
    print 'test_mps'
    N=3
    ms,labels=[],[]
    for i in xrange(N):
        labels.append(('END' if i==0 else 'B%s'%(i-1),'S%s'%i,'END' if i==N-1 else 'B%s'%i))
        if i==0:
            ms.append(array([[[1,0],[0,1]]]))
        elif i==N-1:
            ms.append(array([[[0],[1]],[[1],[0]]]))
        else:
            ms.append(array([[[1,0],[0,1]],[[1,0],[0,1]]]))
    a=MPS(ms,labels)
    #print 'a:\n%s'%a
    print 'a.state: %s'%a.state
    print 'a.norm: %s'%a.norm
    print '-------------------'

    b=deepcopy(a)
    b.canonicalization(cut=b.nsite)
    #print 'b:\n%s'%b
    print 'b.state:%s'%b.state
    print 'b.norm: %s'%b.norm
    print 'b.is_canonical:%s'%(b.is_canonical())
    print '-------------------'

    c=deepcopy(a)
    c.canonicalization(cut=0)
    #print 'c:\n%s'%c
    print 'c.state:%s'%c.state
    print 'c.norm: %s'%c.norm
    print 'c.is_canonical:%s'%(c.is_canonical())
    print

def test_vidal():
    print 'test_vidal'
    #d=c.to_vidal()
    #print 'd:\n%s\n'%d
    #print 'd.state:%s'%d.state
    #print '-------------------'

    #e=d.to_mixed(cut=N)
    #print 'e:\n%s\n'%e
    #print 'e.state:%s'%e.state
    #print '-------------------'
    
    #f=d.to_mixed(cut=0)
    #print 'f:\n%s\n'%f
    #print 'f.state:%s'%f.state
    #print '-------------------'
    print 
