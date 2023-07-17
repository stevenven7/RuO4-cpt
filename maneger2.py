from HamiltonianPy import *
from source import *
from collections import OrderedDict
import numpy as np
import mkl

# tbatasks

def tbatasks(name,parameters,lattice,terms,jobs=()):
    import HamiltonianPy.FreeSystem as TBA
    tba=tbaconstruct(name,parameters,lattice,terms)
    if 'DOS' in jobs:
        tba.register(DOS(name='DOS',ne=401,BZ=KSpace(reciprocals=lattice.reciprocals,nk=100),emin=-5.0,emax=8.0,eta=0.05,savedata=False,plot=True,returndata=False,run=TBA.TBADOS))
    tba.summary()

# vcatasks
def vcatasks(name,parameters,basis,cell,lattice,terms,weiss=(),baths=(),jobs=()):
    import HamiltonianPy.VCA as VCA
    vca=vcaconstruct(name,parameters,[basis],cell,lattice,terms,weiss=())
    if 'EB' in jobs:
        vca.register(VCA.EB(name='EB',path=square_gxm(reciprocals=cell.reciprocals,nk=100),mu=0.0,emin=-5.0,emax=15.0,eta=0.05,ne=401,savedata=False, plot=True,returndata=False,run=VCA.VCAEB))
    if 'DOS' in jobs:
        vca.register(
            DOS(name='DOS', ne=401, BZ=KSpace(reciprocals=cell.reciprocals, nk=100),mu=0.0, emin=-5.0, emax=15.0,eta=0.05,savedata=False, plot=True,returndata=False,run=VCA.VCADOS))
    vca.summary()


if __name__=='__main__':
    mkl.set_num_threads(1)
    Engine.DEBUG=True

    # parameters
    parameters=OrderedDict()
    parameters['E0'] = 1.2
    parameters['E1']=1.3
    parameters['E2']=2.24

    parameters['t0'] = 0.54
    parameters['t1'] = 0.62
    parameters['t2'] = 0.62
    parameters['t12'] = 0

    parameters['U0'] = 1.05
    parameters['U1'] = 1.05
    parameters['U2'] = 1.07
    parameters['U01'] = [0.0,1.03,0.0,0.0]
    parameters['U02'] = [0.0,1.05,0.0,0.0]
    parameters['U12'] = [0.0,1.05,0.0,0.0]

    # tba
    tbatasks(name,parameters,S1('1P-1P',nnb),[E0,E1,E2,t0,t1,t2,t12],jobs=['DOS'])

    # vca
    vcatasks(name, parameters, FBasis(24, 16, 0.0), S1('1P-1P', nnb), S1('2P-2P', nnb),
             [E0, E1, E2, t0, t1, t2, t12, U0, U1, U2, U01, U02, U12], jobs=['EB'])
    vcatasks(name,parameters,FBasis(24,16,0.0),S1('1P-1P',nnb),S1('2P-2P',nnb),[E0,E1,E2,t0,t1,t2,t12,U0,U1,U2,U01,U02,U12],jobs=['DOS'])
