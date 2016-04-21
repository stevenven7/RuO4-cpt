'''
Fermionic terms, including:
1) classes: Quadratic, QuadracticList, Hubbard, HubbardList
2) functions: Hopping, Onsite, Pairing
'''

__all__=['Quadratic','QuadraticList','Hopping','Onsite','Pairing','Hubbard','HubbardList']

from numpy import *
from ..ConstantPy import *
from ..TermPy import *
from ..DegreeOfFreedomPy import *
from ..OperatorPy import * 
from DegreeOfFreedomPy import *
from OperatorPy import * 

class Quadratic(Term):
    '''
    This class provides a complete and unified description for fermionic quadratic terms, i.e. hopping terms, onsite terms and pairing terms.
    Attributes:
        neighbour: integer
            The order of neighbour of this quadratic term.
        indexpackages: IndexPackageList or function which returns IndexPackageList
            The indexpackages of the quadratic term.
            When it is a function, it can return bond dependent indexpackages as needed.
        amplitude: function which returns float or complex
            This function can return bond dependent and index dependent coefficient as needed.
    Note: The final coefficient comes from three parts, the value of itself, the value of the indexpacakge, and the value amplitude returns.
    '''
    
    def __init__(self,id,mode,value,neighbour=0,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
        '''
        Constructor.
        Parameters:
            id: string
                The specific id of the term.
            mode: string
                The type of the term.
            value: float or complex
                The overall coefficient of the term.
            neighbour: integer, optional
                The order of neighbour of the term.
            atoms,orbitals,spins: list of integer, optional
                The atom, orbital and spin index used to construct a wanted instance of IndexPackage.
                When the parameter indexpackages is a function, these parameters will be omitted.
            indexpackages: IndexPackageList or function
                1) IndexPackageList:
                    It will be multiplied by an instance of IndexPackage constructed from atoms, orbitals and spins as the final indexpackages. 
                2) function:
                    It must return an instance of IndexPackageList and take an instance of Bond as its only argument.
            amplitude: function
                It must return a float or complex and take an instance of Bond as its only argument.
            modulate: function
                It must return a float or complex and its arguments are unlimited.
        '''
        super(Quadratic,self).__init__(id,mode,value,modulate)
        self.neighbour=neighbour
        if indexpackages is not None:
            if isinstance(indexpackages,IndexPackageList):
                self.indexpackages=IndexPackage(1,atoms=atoms,orbitals=orbitals,spins=spins)*indexpackages
            elif callable(indexpackages):
                self.indexpackages=indexpackages
            else:
                raise ValueError('Quadratic init error: the input indexpackages should be an instance of IndexPackageList or a function.')
        else:
            self.indexpackages=IndexPackageList(IndexPackage(1,atoms=atoms,orbitals=orbitals,spins=spins))
        self.amplitude=amplitude

    def __repr__(self):
        '''
        Convert an instance to string.
        '''
        result=[]
        result.append('id=%s'%self.id)
        result.append('mode=%s'%self.mode)
        result.append('value=%s'%self.value)
        result.append('neighbour=%s'%self.neighbour)
        result.append('indexpackages=%s'%self.indexpackages)
        if self.amplitude is not None:
            result.append('amplitude=%s'%self.amplitude)
        if hasattr(self,'modulate'):
            result.append('modulate=%s'%self.modulate)
        return 'Quadratic('+', '.join(result)+')'

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a Quadratic instance with a Quadratic/QuadraticList instance.
        '''
        result=QuadraticList()
        result.append(self)
        if isinstance(other,Quadratic):
            result.append(other)
        elif isinstance(other,QuadraticList):
            result.extend(other)
        else:
            raise ValueError('Quadratic "+" error: the other parameter must be an instance of Quadratic or QuadraticList.')
        return result

    def __pos__(self):
        '''
        Overloaded operator(+), i.e. +self.
        '''
        result=QuadraticList()
        result.append(self)
        return result

    def mesh(self,bond,config,dtype=complex128):
        '''
        This method returns the mesh of a quadratic term defined on a bond.
        Parameters:
            bond: Bond
                The bond on which the quadratic term is defined.
            config: Configuration
                The configuration of degrees of freedom.
            dtype: complex128,complex64,float64,float32, optional
                The data type of the returned mesh.
        Returns: 2D ndarray
            The mesh.
        '''
        edgr=config[bond.epoint.pid]
        sdgr=config[bond.spoint.pid]
        n1=edgr.norbital*edgr.nspin*edgr.nnambu
        n2=sdgr.norbital*sdgr.nspin*sdgr.nnambu
        result=zeros((n1,n2),dtype=dtype)
        if self.neighbour==bond.neighbour:
            value=self.value*(1 if self.amplitude==None else self.amplitude(bond))
            if callable(self.indexpackages):
                buff=self.indexpackages(bond)
            else:
                buff=self.indexpackages
            for obj in buff:
                pv=value*obj.value
                eatom=edgr.atom
                satom=sdgr.atom
                if hasattr(obj,'atoms'):
                    eatom=obj.atoms[0]
                    satom=obj.atoms[1]
                if eatom==edgr.atom and satom==sdgr.atom:
                    enambu=CREATION if self.mode=='pr' else ANNIHILATION
                    snambu=ANNIHILATION
                    if hasattr(obj,'spins'):
                        if hasattr(obj,'orbitals'):
                            i=edgr.seq_state(FID(obj.orbitals[0],obj.spins[0],enambu))
                            j=sdgr.seq_state(FID(obj.orbitals[1],obj.spins[1],snambu))
                            result[i,j]+=pv
                        elif edgr.norbital==sdgr.norbital:
                            for k in xrange(edgr.norbital):
                                i=edgr.seq_state(FID(k,obj.spins[0],enambu))
                                j=sdgr.seq_state(FID(k,obj.spins[1],snambu))
                                result[i,j]+=pv
                        else:
                            raise ValueError("Quadratic mesh error: the norbital on the epoint and the spoint of the input bond should be equal.")
                    elif edgr.nspin==sdgr.nspin:
                        if hasattr(obj,'orbitals'):
                            for k in xrange(edgr.nspin):
                                i=edgr.seq_state(FID(obj.orbitals[0],k,enambu))
                                j=sdgr.seq_state(FID(obj.orbitals[1],k,snambu))
                                result[i,j]+=pv
                        elif n1==n2:
                            ns=edgr.norbital*edgr.nspin
                            if self.mode=='pr':
                                for k in xrange(ns):
                                    result[k,k+ns]+=pv
                            else:
                                for k in xrange(ns):
                                    result[k,k]+=pv
                        else:
                            raise ValueError("Quadratic mesh error: the norbital of epoint and spoint of the input bond should be equal.")
                    else:
                        raise ValueError("Quadratic mesh error: the nspin of epoint and spoint of the input bond should be equal.")
            if self.neighbour==0:
                for i in xrange(n1):
                    result[i,i]/=2.0
                    for j in xrange(i):
                        if abs(result[i,j]-conjugate(result[j,i]))<RZERO: result[i,j]=0
        return result

def Hopping(id,value,neighbour=1,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
    '''
    A specified function to construct a hopping term.
    '''
    return Quadratic(id,'hp',value,neighbour,atoms,orbitals,spins,indexpackages,amplitude,modulate)

def Onsite(id,value,neighbour=0,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
    '''
    A specified function to construct an onsite term.
    '''
    return Quadratic(id,'st',value,neighbour,atoms,orbitals,spins,indexpackages,amplitude,modulate)

def Pairing(id,value,neighbour=0,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
    '''
    A specified function to construct an pairing term.
    '''
    return Quadratic(id,'pr',value,neighbour,atoms,orbitals,spins,indexpackages,amplitude,modulate)

class QuadraticList(TermList):
    '''
    This class packs several Quadratic instances as a whole for convenience.
    '''

    def mesh(self,bond,config,mask=None,dtype=complex128):
        '''
        This method returns the mesh of all quadratic terms defined on a bond.
        Parameters:
            bond: Bond
                The bond on which the quadratic terms are defined.
            config: Configuration
                The configuration of degrees of freedom.
            mask: callable
                A function used to pick the quadratic terms whose only argument is an instance of Quadratic. 
                If the returned value if True, the masked quadratic term is included.
            dtype: complex128,complex64,float128,float64, optional
                The data type of the returned mesh.
        Returns: 2D ndarray
            The mesh.
        '''
        edgr,sdgr=config[bond.epoint.pid],config[bond.spoint.pid]
        if edgr.nnambu==sdgr.nnambu:
            result=zeros((edgr.norbital*edgr.nspin*edgr.nnambu,sdgr.norbital*sdgr.nspin*sdgr.nnambu),dtype=dtype)
            for obj in self:
                if mask is None or mask(obj):
                    result+=obj.mesh(bond,config,dtype=dtype)
            return result
        else:
            raise ValueError('Quadratic mesh error: the nnambu of epoint and spoint must be equal.')

    def operators(self,bond,table,config,mask=None,dtype=complex128):
        '''
        This method returns all the desired quadratic operators defined on the input bond with non-zero coefficients.
        Parameters:
            bond: Bond
                The bond where the quadratic operators are defined.
            table: Table
                The index-sequence table.
                Only those operators whose indices are in this table will be returned.
            config:
                The configuration of degrees of freedom.
            mask: callable
                A function used to pick the quadratic terms whose only argument is an instance of Quadratic. 
                If the returned value if True, the operators corresponding to the masked quadratic term are included.
            dtype: dtype, optional
                The data type of the coefficient of the returned operators.
        Returns:
            result: OperatorCollection
                The quadratic operators with non-zero coefficients.
        Note: only one "half" of the operators are returned, which means
                1) the Hermitian conjugate of non-Hermitian operators is not included;
                2) the coefficient of the self-Hermitian operators is divided by a factor 2;
                3) the BdG case, only the electron part of the hopping terms and onsite terms are contained, and for the electron part, Rule 1) and 2) also applies.
        '''
        result=to_operators(self.mesh(bond,config,mask,dtype=dtype),bond,table,config)
        if bond.neighbour!=0:
            result+=to_operators(self.mesh(bond.reversed,config,mask=lambda quadratic: True if (quadratic.mode=='pr' and (True if mask is None else mask(quadratic))) else False,dtype=dtype),bond.reversed,table,config)
        return result

def to_operators(mesh,bond,table,config):
    result=OperatorCollection()
    indices=argwhere(abs(mesh)>RZERO)
    for (i,j) in indices:
        eindex=Index(bond.epoint.pid,config[bond.epoint.pid].state_index(i))
        sindex=Index(bond.spoint.pid,config[bond.spoint.pid].state_index(j))
        if eindex in table and sindex in table:
            result+=F_Quadratic(mesh[i,j],indices=(eindex.replace(nambu=1-eindex.nambu),sindex),rcoords=[bond.rcoord],icoords=[bond.icoord],seqs=(table[eindex],table[sindex]))
    return result

class Hubbard(Term):
    '''
    This class provides a complete description of single orbital and multi orbital Hubbard interactions.
    Attributes:
        value: float or 1D array-like, inherited from Term
            float: single-orbital Hubbard interaction.
            1D array-like: multi-orbital Hubbard interaction and the length must be 3.
                value[0]: intra-orbital interaction 
                value[1]: inter-orbital interaction
                value[2]: spin-flip and pair-hopping interaction
        atom: integer 
            The atom index of the point where the Hubbard interactions are defined.
    '''

    def __init__(self,id,value,atom=None,modulate=None):
        '''
        Constructor.
        '''
        super(Hubbard,self).__init__('hb',id,value,modulate)
        if atom is not None: self.atom=atom

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=[]
        result.append('id=%s'%self.id)
        if hasattr(self,'atom'): 
            result.append('atom=%s'%self.atom)
        result.append('value=%s'%self.value)
        return 'Hubbard('+', '.join(result)+')'

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a Hubbard instance with a Hubbard/HubbardList instance.
        '''
        result=HubbardList()
        result.append(self)
        if isinstance(other,Hubbard):
            result.append(other)
        elif isinstance(other,HubbardList):
            result.extend(other)
        else:
            raise ValueError('Hubbard "+" error: the other parameter must be an instance of Hubbard or HubbardList')
        return result

    def __pos__(self):
        '''
        Overloaded operator(+), i.e. +self.
        '''
        result=HubbardList()
        result.append(self)
        return result

    def mesh(self,bond,config,dtype=float64):
        '''
        This method returns the mesh of Hubbard terms.
        Parameters:
            bond: Bond
                The bond on which the Hubbard term is defined.
            config: Configuration
                The configuration of degrees of freedom.
            dtype: complex128,complex64,float128,float64, optional
                The data type of the returned mesh.
        Returns: 4D ndarray
            The mesh.
        '''
        dgr=config[bond.epoint.pid]
        ndim=dgr.norbital*dgr.nspin
        result=zeros((ndim,ndim,ndim,ndim),dtype=dtype)
        if hasattr(self,'atom'):
            atom=self.atom
        else:
            atom=dgr.atom
        if atom==dgr.atom:
            try:
                nv=len(self.value)
            except TypeError:
                nv=1
            if nv>=1:
                for h in xrange(dgr.norbital):
                    i=dgr.seq_state(FID(h,1,ANNIHILATION))
                    j=dgr.seq_state(FID(h,0,ANNIHILATION))
                    k=j
                    l=i
                    result[i,j,k,l]=self.value/2 if nv==1 else self.value[0]/2
            if nv==3:
                for h in xrange(dgr.norbital):
                    for g in xrange(dgr.norbital):
                      if g!=h:
                        i=dgr.seq_state(FID(g,1,ANNIHILATION))
                        j=dgr.seq_state(FID(h,0,ANNIHILATION))
                        k=j
                        l=i
                        result[i,j,k,l]=self.value[1]/2
                for h in xrange(dgr.norbital):
                    for g in xrange(h):
                        for f in xrange(2):
                            i=dgr.seq_state(FID(g,f,ANNIHILATION))
                            j=dgr.seq_state(FID(h,f,ANNIHILATION))
                            k=j
                            l=i
                            result[i,j,k,l]=(self.value[1]-self.value[2])/2
                for h in xrange(dgr.norbital):
                    for g in xrange(h):
                        i=dgr.seq_state(FID(g,1,ANNIHILATION))
                        j=dgr.seq_state(FID(h,0,ANNIHILATION))
                        k=dgr.seq_state(FID(g,0,ANNIHILATION))
                        l=dgr.seq_state(FID(h,1,ANNIHILATION))
                        result[i,j,k,l]=self.value[2]
                for h in xrange(dgr.norbital):
                    for g in xrange(h):
                        i=dgr.seq_state(FID(g,1,ANNIHILATION))
                        j=dgr.seq_state(FID(g,0,ANNIHILATION))
                        k=dgr.seq_state(FID(h,0,ANNIHILATION))
                        l=dgr.seq_state(FID(h,1,ANNIHILATION))
                        result[i,j,k,l]=self.value[2]
        return result

class HubbardList(list):
    '''
    This class pack several Hubbard instances as a whole for convenience.
    '''

    def mesh(self,bond,config,dtype=float64):
        '''
        This method returns the mesh of all Hubbard terms defined on a bond.
        Parameters:
            bond: Bond
                The bond on which the Hubbard term is defined.
            config: Configuration
                The configuration of degrees of freedom.
            dtype: complex128,complex64,float128,float64, optional
                The data type of the returned mesh.
        Returns: 4D ndarray
            The mesh.
        '''
        if bond.neighbour==0:
            edgr,sdgr=config[bond.epoint.pid],config[bond.spoint.pid]
            if edgr.nspin==2 and sdgr.nspin==2:
                ndim=edgr.norbital*edgr.nspin
                result=zeros((ndim,ndim,ndim,ndim),dtype=dtype)
                for obj in self:
                    result+=obj.mesh(bond,config,dtype=dtype)
                return result
            else:
                raise ValueError('HubbardList mesh error: the input bond must be onsite ones nspin=2.')
        else:
            return 0

    def operators(self,bond,table,config,dtype=float64):
        '''
        This method returns all the Hubbard operators defined on the input bond with non-zero coefficients.
        Parameters:
            bond: Bond
                The bond on which the Hubbard terms are defined.
            table: Table
                The index-sequence table.
                Since Hubbard terms are quartic, it never uses the Nambu space.
            config: Configuration
                The configuration of degrees of freedom.
            dtype: dtype,optional
                The data type of the coefficient of the returned operators.
        Returns:
            result: OperatorCollection
                All the Hubbard operators with non-zero coeffcients.
        Note: only one "half" of the operators are returned, which means
                1) the Hermitian conjugate of non-Hermitian operators is not included;
                2) the coefficient of the self-Hermitian operators is divided by a factor 2.
        '''
        result=OperatorCollection()
        dgr=config[bond.epoint.pid]
        mesh=self.mesh(bond,config,dtype=dtype)
        indices=argwhere(abs(mesh)>RZERO)
        for (i,j,k,l) in indices:
            index1=Index(bond.epoint.pid,dgr.state_index(i))
            index2=Index(bond.epoint.pid,dgr.state_index(j))
            index3=Index(bond.epoint.pid,dgr.state_index(k))
            index4=Index(bond.epoint.pid,dgr.state_index(l))
            result+=F_Hubbard(mesh[i,j,k,l],indices=(index1.replace(nambu=CREATION),index2.replace(nambu=CREATION),index3,index4),rcoords=[bond.epoint.rcoord],icoords=[bond.epoint.icoord],seqs=(table[index1],table[index2],table[index3],table[index4]))
        return result
