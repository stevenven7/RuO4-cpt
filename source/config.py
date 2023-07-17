from HamiltonianPy import *
import numpy as np

__all__=['name','nnb','parametermap','idfmap','E0','E1','E2','t0','t1','t2','t12','U0','U1','U2','U01','U02','U12','S1']


# The configs of the model
name='yoodee'
nnb=2

# parametermap
parametermap=None

# idfmap
idfmap=lambda pid: Fock(atom=pid.site%1.5,norbital=3,nspin=2,nnambu=1)

#amplitude(bond)
def amplitudeE0(bond,config,index1=None,index2=None):
    if index1.orbital==index2.orbital==0:
        result=1.0
    else:
        result=0.0
    return result

def amplitudeE1(bond,config,index1=None,index2=None):
    if index1.orbital==index2.orbital==1:
        result=1.0
    else:
        result=0.0
    return result

def amplitudeE2(bond,config,index1=None,index2=None):
    if index1.orbital==index2.orbital==2:
        result=1.0
    else:
        result=0.0
    return result

def amplitudet0(bond,config,index1=None,index2=None):
    if index1.orbital==index2.orbital==0:
        result=1.0
    else:
        result=0.0
    return result

def amplitudet1(bond,config,index1=None,index2=None):
    if index1.orbital==index2.orbital==1:
        result=1.0
    else:
        result=0.0
    return result

def amplitudet2(bond,config,index1=None,index2=None):
    if index1.orbital==index2.orbital==2:
        result=1.0
    else:
        result=0.0
    return result

def amplitudet12(bond,config,index1=None,index2=None):
    if index1.orbital==1 and index2.orbital==2:
        result=1.0
    elif index1.orbital==2 and index2.orbital==1:
        result = 1.0
    else:
        result=0.0
    return result

def amplitudeU0(bond,config,index1=None,index2=None):
    if index1==None:
        result=1.0
    elif index1.orbital==index2.orbital==0:
        result=1.0
    else:
        result=0.0
    return result

def amplitudeU1(bond,config,index1=None,index2=None):
    if index1==None:
        result=1.0
    elif index1.orbital==index2.orbital==1:
        result=1.0
    else:
        result=0.0
    return result

def amplitudeU2(bond,config,index1=None,index2=None):
    if index1==None:
        result=1.0
    elif index1.orbital==index2.orbital==2:
        result=1.0
    else:
        result=0.0
    return result

def amplitudeU01(bond,config,index1=None,index2=None):
    if index1==None:
        result=1.0
    elif index1.orbital==0 and index2.orbital==1:
        result=1.0
    elif index1.orbital==1 and index2.orbital ==0:
        result = 1.0
    else:
        result=0.0
    return result

def amplitudeU02(bond,config,index1=None,index2=None):
    if index1==None:
        result=1.0
    elif index1.orbital== 0 and index2.orbital== 2:
        result=1.0
    elif index1.orbital== 2 and index2.orbital == 0:
        result = 1.0
    else:
        result=0.0
    return result

def amplitudeU12(bond,config,index1=None,index2=None):
    if index1==None:
        result=1.0
    elif index1.orbital== 1 and index2.orbital== 2:
        result = 1.0
    elif index1.orbital == 2 and index2.orbital == 1:
        result = 1.0
    else:
        result=0.0
    return result

# terms

E0=lambda **parameters: Onsite('E0',parameters['E0'],amplitude=amplitudeE0)
E1=lambda **parameters: Onsite('E1',parameters['E1'],amplitude=amplitudeE1)
E2=lambda **parameters: Onsite('E2',parameters['E2'],amplitude=amplitudeE2)

t0=lambda **parameters: Hopping('t0',parameters['t1'],neighbour=1,amplitude=amplitudet0)
t1=lambda **parameters: Hopping('t1',parameters['t1'],neighbour=1,amplitude=amplitudet1)
t2=lambda **parameters: Hopping('t2',parameters['t1'],neighbour=1,amplitude=amplitudet2)
t12=lambda **parameters: Hopping('t12',parameters['t1'],neighbour=1,amplitude=amplitudet12)

U0=lambda **parameters: Hubbard('U0',parameters['U0'],amplitude=amplitudeU0)
U1=lambda **parameters: Hubbard('U1',parameters['U1'],amplitude=amplitudeU1)
U2=lambda **parameters: Hubbard('U2',parameters['U2'],amplitude=amplitudeU2)
U01=lambda **parameters: Hubbard('U01',parameters['U01'],amplitude=amplitudeU01)
U02=lambda **parameters: Hubbard('U02',parameters['U02'],amplitude=amplitudeU02)
U12=lambda **parameters: Hubbard('U12',parameters['U12'],amplitude=amplitudeU12)

# cluster
S1 = Square('S1')
