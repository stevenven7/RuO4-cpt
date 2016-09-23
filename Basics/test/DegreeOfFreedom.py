'''
Degree of freedom test.
'''

__all__=['test_deg_fre']

from HamiltonianPy.Basics.Geometry import *
from HamiltonianPy.Basics.Tree import *
from HamiltonianPy.Basics.DegreeOfFreedom import *
from HamiltonianPy.Basics.FermionicPackage.DegreeOfFreedom import *
from HamiltonianPy.Basics.SpinPackage.DegreeOfFreedom import *

def test_deg_fre():
    test_table()
    test_deg_fre_tree()

def test_table():
    print 'test_table'
    a=Table(['i1','i2'])
    b=Table(['i3','i4'])
    c=Table.union([a,b],key=lambda key: key[1])
    print 'a: %s'%a
    print 'b: %s'%b
    print 'c=union(a,b): %s'%c
    print 'c.reverse_table: %s'%c.reversed_table
    print 'c["i4"]: %s'%c['i4']
    print 'c.subset: %s'%c.subset(mask=lambda key: True if key!='i1' else False)
    print

def test_deg_fre_tree():
    print 'test_deg_fre_tree'
    config=Configuration(priority=DEFAULT_FERMIONIC_PRIORITY)
    for site in xrange(4):
        config[PID(scope=1,site=site)]=Fermi(norbital=2,nspin=2,nnambu=1)
        config[PID(scope=2,site=site)]=Fermi(norbital=2,nspin=2,nnambu=1)
    tree=DegFreTree(config,layers=DEFAULT_FERMI_LAYERS)
    for layer in DEFAULT_FERMI_LAYERS:
        print 'layer: %s'%layer
        for i,index in enumerate(tree.indices(layer)):
            print i,index,tree[index]
        print
    config=Configuration(priority=DEFAULT_SPIN_PRIORITY)
    for site in xrange(4):
        config[PID(scope=1,site=site)]=Spin(S=0.5)
        config[PID(scope=2,site=site)]=Spin(S=0.5)
    tree=DegFreTree(config,layers=DEFAULT_SPIN_LAYERS)
    for layer in DEFAULT_SPIN_LAYERS:
        print 'layer: %s'%layer
        for i,index in enumerate(tree.indices(layer)):
            print i,index,tree[index]
        print
    print
