'''
Degrees of freedom in a lattice, including:
1) classes: Table, Index, Internal, Configuration, DegFreTree, IndexPack, IndexPackList
'''

__all__=['Table','Index','Internal','Configuration','DegFreTree','IndexPack','IndexPackList']

import numpy as np
from Tree import Tree
from Geometry import PID
from collections import OrderedDict

class Table(dict):
    '''
    This class provides the methods to get an index from its sequence number or vice versa.
    '''
    def __init__(self,indices=[],key=None):
        '''
        Constructor.
        Parameters:
            indices: list of any hashable object
                The indices that need to be mapped to sequences.
            key: function, optional
                The function used to sort the indices.
            NOTE: The final order of the index in indices will be used as its sequence number.
        '''
        for i,v in enumerate(indices if key is None else sorted(indices,key=key)):
            self[v]=i

    @staticmethod
    def union(tables,key=None):
        '''
        This function returns the union of index-sequence tables.
        Parameters:
            tables: list of Table
                The tables to be unioned.
            key: callable, optional
                The function used to compare different indices in tables.
                When it is None, the sequence of an index will be naturally ordered by the its sequence in the input tables.
        Returns: Table
            The union of the input tables.
        '''
        result=Table()
        if key is None:
            sum=0
            for table in tables:
                if isinstance(table,Table):
                    count=0
                    for k,v in table.iteritems():
                        result[k]=v+sum
                        count+=1
                    sum+=count
        else:
            for table in tables:
                result.update(table)
            for i,k in enumerate(sorted(result,key=key)):
                result[k]=i
        return result

    def subset(self,mask):
        '''
        This function returns a certain subset of an index-sequence table according to the mask function.
        Parameters:
            mask: callable
                A certain subset of table is extracted according to the return value of this function on the index in the table.
                When the return value is True, the index will be included and the sequence is naturally determined by its order in the mother table.
        Returns:
            The subset table.
        '''
        result=Table()
        for k,v in self.iteritems():
            if mask(k):
                result[k]=v
        buff={}
        for i,k in enumerate(sorted([key for key in result.keys()],key=result.get)):
            buff[k]=i
        result.update(buff)
        return result

    @property
    def reversed_table(self):
        '''
        This function returns the sequence-index table for a reversed lookup.
        Returns: Table
            The reversed table whose key is the sequence and value the index.
        '''
        result=Table()
        for k,v in self.iteritems():
            result[v]=k
        return result

class Index(tuple):
    '''
    This class provides an index for a microscopic degree of freedom, including the spatial part and interanl part.
    '''
    def __new__(cls,pid,iid):
        '''
        Constructor.
        Parameters:
            pid: PID
                The point index, i.e. the spatial part in a lattice of the index
            iid: namedtuple
                The internal index, i.e. the internal part in a point of the index.
        '''
        self=super(Index,cls).__new__(cls,pid+iid)
        self.__dict__=OrderedDict()
        self.__dict__.update(pid._asdict())
        self.__dict__.update(iid._asdict())
        return self

    def __copy__(self):
        '''
        Copy.
        '''
        return self.replace(**self.__dict__)

    def __deepcopy__(self,memo):
        '''
        Deep copy.
        '''
        return self.replace(**self.__dict__)

    @property
    def pid(self):
        '''
        The pid part of the index.
        '''
        return PID(**{key:getattr(self,key) for key in PID._fields})

    def iid(self,cls):
        '''
        The iid part of the index.
        '''
        return cls(**{key:value for key,value in self.__dict__.items() if key not in PID._fields})

    def __repr__(self):
        '''
        Convert an instance to string.
        '''
        return ''.join(['Index(','=%r, '.join(self.__dict__.keys()),'=%r)'])%self

    def replace(self,**karg):
        '''
        Return a new Index object with specified fields replaced with new values.
        '''
        result=tuple.__new__(Index,map(karg.pop,self.__dict__.keys(),self))
        if karg:
            raise ValueError('Index replace error: it got unexpected field names: %r'%karg.keys())
        result.__dict__=OrderedDict()
        for key,value in zip(self.__dict__.keys(),result):
            result.__dict__[key]=value
        return result

    def to_tuple(self,priority):
        '''
        Convert an instance to tuple according to the parameter priority.
        Parameters:
            priority: list of string
                Every element of this list should correspond to a name of an attribute of self.
                The elements should have no duplicates and its length should be equal to the number of self's attributes.
        Returns: tuple
            The elements of the returned tuple are equal to the attributes of self. 
            Their orders are determined by the orders they appear in priority.
        '''
        if len(set(priority))<len(priority):
            raise ValueError('Index to_tuple error: the priority has duplicates.')
        if len(priority)!=len(self.__dict__):
            raise ValueError("Index to_tuple error: the priority doesn't cover all the attributes.")
        return tuple(map(self.__dict__.get,priority))

class Internal(object):
    '''
    This class is the base class for all internal degrees of freedom in a single point.
    '''

    def ndegfre(self,mask=None):
        '''
        Return the number of the interanl degrees of freedom modified by mask.
        Parameters:
            mask: list of string, optional
                Only the indices in mask can be varied in the counting of the number of the degrees of freedom.
                When None, all the allowed indices can be varied and thus the total number of the interanl degrees of freedom is returned.
        Returns: number
            The requested number of the interanl degrees of freedom.
        '''
        raise NotImplementedError()

    def indices(self,pid,**karg):
        '''
        Return a list of all the allowed indices within this internal degrees of freedom combined with an extra spatial part.
        Parameters:
            pid: PID
                The extra spatial part of the indices.
        Returns: list of Index
            The allowed indices.
        '''
        raise NotImplementedError("%s indices error: it is not implemented."%self.__class__.__name__)

class Configuration(dict):
    '''
    Configuration of the degrees of freedom in a lattice.
    Attributes:
        priority: list of string
            The sequence priority of the allowed indices that can be defined on a lattice.
    '''

    def __init__(self,dict=None,priority=None):
        '''
        Constructor.
        Parameters:
            dict: dict with key PID and value Internal respectively
                The key is the pid of the lattice point and the value is the internal degrees of freedom of that point.
            priority: list of string
                The sequence priority of the allowed indices that can be defined on the lattice.
        '''
        if dict is not None:
            self.update(dict)
        self.priority=priority

    def table(self,**karg):
        '''
        Return a Table instance that contains all the allowed indices which can be defined on a lattice.
        '''
        return Table([index for key,value in self.items() for index in value.indices(key,**karg)],key=lambda index: index.to_tuple(priority=self.priority))

    def __add__(self,other):
        '''
        Add two configurations.
        '''
        assert self.priority==other.priority
        result=Configuration(priority=self.priority)
        result.update(self)
        result.update(other)
        return result

class DegFreTree(Tree):
    '''
    The tree of the layered degrees of freedom.
    For each (node,data) pair of the tree,
        node: Index
            The masked index which can represent a couple of indices.
        data: integer
            The number of degrees of freedom that the index represents.
    Attributes:
        config: Configuration
            The configuration of the interanl degrees of freedom in a lattice.
        layers: list of string
            The tag of each layer of indices.
    '''

    def __init__(self,config,layers=None):
        '''
        Constructor.
        Parameters:
            config: Configuration
                The configuration of the interanl degrees of freedom in a lattice.
            layers: list of string
                The tag of each layer of indices.
        '''
        assert set(list(xrange(len(PID._fields))))==set([layers.index(key) for key in PID._fields])
        Tree.__init__(self,root=tuple([None]*len(layers)),data=None)
        indices=config.table()
        for i,layer in enumerate(layers):
            leaves=set([index.replace(**{key:None for key in layers[i+1:]}) for index in indices])
            for leaf in leaves:
                self.add_leaf(parent=leaf.replace(**{layer:None}),leaf=leaf,data=None)
        for i,indices in enumerate(self.indices(*list(xrange(len(layers)-1,-1,-1)))):
            for index in indices:
                if i>len(layers)-len(PID._fields):
                    self[index]=np.product([self[child] for child in self.children(index)])
                else:
                    self[index]=config[index.pid].ndegfre(mask=[key for key in layers[len(PID._fields):] if getattr(index,key) is None])
        self.config=config
        self.layers=layers

    def indices(self,*layers):
        '''
        The indices in some layers.
        Parameters:
            layers: list of string/integer
                The tags/numbers of the requested layers of indices.
        Returns: list of list of Index when len(layers)==1 else list of Index
            The indices in the requested layers.
        '''
        result,buff=[None]*len(layers),[self.root]
        layers=[(self.layers.index(layer) if isinstance(layer,str) else layer) for layer in layers]
        for i in xrange(max(layers)+1):
            temp=[]
            for key in buff[:]:
                temp.extend(self.children(key))
            buff=temp
            if i in layers:
                if len(layers)==1:
                    result=buff
                else:
                    result[layers.index(i)]=buff
        return result

class IndexPack(object):
    '''
    This class packs several degrees of freedom as a whole for convenience.
    Attributes:
        value: float64 or complex128
            The overall coefficient of the IndexPack.
    '''

    def __init__(self,value):
        '''
        Constructor.
        Parameters:
            value: float64/complex128
                The overall coefficient of the IndexPack.
        '''
        self.value=value

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of an IndexPack instance with an IndexPack/IndexPackList instance.
        '''
        result=IndexPackList()
        result.append(self)
        if issubclass(other.__class__,IndexPack):
            result.append(other)
        elif isinstance(other,IndexPackList):
            result.extend(other)
        else:
            raise ValueError("IndexPack '+' error: the 'other' parameter must be of class IndexPack or IndexPackList.")
        return result

    def __rmul__(self,other):
        '''
        Overloaded operator(*).
        '''
        return self.__mul__(other)

class IndexPackList(list):
    '''
    This class packs several IndexPack as a whole for convenience.
    '''

    def __init__(self,*arg):
        for buff in arg:
            if issubclass(buff.__class__,IndexPack):
                self.append(buff)
            else:
                raise ValueError("IndexPackList init error: the input parameters must be of IndexPack's subclasses.")

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        return 'IndexPackList('+', '.join([str(obj) for obj in self])

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of an IndexPackList instance with an IndexPack/IndexPackList instance.
        '''
        result=IndexPackList(*self)
        if isinstance(other,IndexPack):
            result.append(other)
        elif isinstance(other,IndexPackList):
            result.extend(other)
        else:
            raise ValueError("IndexPackList '+' error: the 'other' parameter must be of class IndexPack or IndexPackList.")
        return result

    def __radd__(self,other):
        '''
        Overloaded operator(+).
        '''
        return self.__add__(other)

    def __mul__(self,other):
        '''
        Overloaded operator(*).
        '''
        result=IndexPackList()
        for buff in self:
            temp=buff*other
            if isinstance(temp,IndexPackList):
                result.extend(temp)
            elif issubclass(temp.__class__,IndexPack):
                result.append(temp)
            else:
                raise ValueError("IndexPackList *' error: the element(%s) in self multiplied by other is not of IndexPack/IndexPackList."%(buff))
        return result

    def __rmul__(self,other):
        '''
        Overloaded operator(*).
        '''
        return self.__mul__(other)
