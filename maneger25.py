from HamiltonianPy import *
from source import *
from collections import OrderedDict
import numpy as np
import mkl

# tbatasks

def tbatasks(name,parameters,lattice,terms,jobs=()):
    import HamiltonianPy.FreeSystem as TBA
    tba=tbaconstruct(name,parameters,lattice,terms)
    if 'EB' in jobs:
        if len(lattice.vectors)==2:
            tba.register(EB(name='EB',path=square_gxm(reciprocals=lattice.reciprocals,nk=100),run=TBA.TBAEB))
    if 'DOS' in jobs:
        tba.register(DOS(name='DOS',ne=401,BZ=KSpace(reciprocals=lattice.reciprocals,nk=100),emin=-5.0,emax=5.0,eta=0.05,savedata=False,plot=True,returndata=False,run=TBA.TBADOS))
    tba.summary()


# edtasks
def edtasks(name,parameters,basis,lattice,terms,jobs=()):
    import HamiltonianPy.ED as ED
    ed=edconstruct(name,parameters,[basis],lattice,terms)
    if 'EL' in jobs: ed.register(ED.EL(name='EL',path=BaseSpace(['U',np.linspace(0,30.0,301)]),ns=1,nder=2,run=ED.EDEL))
    ed.summary()

# vcatasks
def vcatasks(name,parameters,basis,cell,lattice,terms,weiss=(),baths=(),jobs=()):
    import HamiltonianPy.VCA as VCA
    vca=vcaconstruct(name,parameters,[basis],cell,lattice,terms,weiss=())
    if 'EB' in jobs:
        vca.register(VCA.EB(name='EB',path=square_gxm(reciprocals=cell.reciprocals,nk=100),mu=0.0,emin=3.0,emax=12.0,eta=0.05,ne=401,savedata=False, plot=True,returndata=False,run=VCA.VCAEB))
    if 'DOS' in jobs:
        vca.register(
            DOS(name='DOS', ne=401, BZ=KSpace(reciprocals=cell.reciprocals, nk=100),mu=0.0, emin=3.0, emax=12.0,eta=0.05,savedata=False, plot=True,returndata=False,run=VCA.VCADOS))
    vca.summary()


if __name__=='__main__':
    mkl.set_num_threads(1)
    Engine.DEBUG=True

    # parameters
    parameters=OrderedDict()
    parameters['E0'] = 1.2
    parameters['E1']=1.3
    parameters['E2']=2.24

    parameters['t0'] = 0.617
    parameters['t1'] = 0.685
    parameters['t2'] = 0.686
    parameters['t12'] = 0

    parameters['U0'] = 2.20
    parameters['U1'] = 1.58
    parameters['U2'] = 1.61
    parameters['U01'] = [0.0,1.86,0.0,0.0]
    parameters['U02'] = [0.0,1.90,0.0,0.0]
    parameters['U12'] = [0.0,1.59,0.0,0.0]

    # tba
    tbatasks(name, parameters, S1('1P-1P', nnb), [E0, E1, E2, t0, t1, t2, t12], jobs=['EB'])
    tbatasks(name,parameters,S1('1P-1P',nnb),[E0,E1,E2,t0,t1,t2,t12],jobs=['DOS'])

    # ed
    #edtasks(name,parameters,FBasis(24,16,0.0),S1('1P-1P', nnb),[E0,E1,E2,t0,t1,t2,t12,U0,U1,U2,U01,U02,U12],jobs=['EL'])

    # vca
    vcatasks(name, parameters, FBasis(24, 16, 0.0), S1('1P-1P', nnb), S1('2P-2P', nnb),
             [E0, E1, E2, t0, t1, t2, t12, U0, U1, U2, U01, U02, U12], jobs=['EB'])
    vcatasks(name,parameters,FBasis(24,16,0.0),S1('1P-1P',nnb),S1('2P-2P',nnb),[E0,E1,E2,t0,t1,t2,t12,U0,U1,U2,U01,U02,U12],jobs=['DOS'])
