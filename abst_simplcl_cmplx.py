## stdlib
import collections

## 3rd party
import numpy as np

# Team : Bill Lee, Raghavendra Padmanabhan and Francisco Vargas
#--------------------------------------------------------------------------------------

class Vertex:
    ''' Representation of a vertex in the ASC class
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.label            = kwargs["label"]
        self.parent           = kwargs.get("parent", None)
        self.level            = kwargs.get("level", None)
        self.connections      = kwargs.get("connections", set()) # child connections

    def __repr__(self):
        '''
        '''
        return "Vertex(label = {}, level = {})".format(self.label, self.level)


class ASCInsertionError(Exception):
    ''' Exception raised for invalid insertions to ASC
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.simplex_inserted = kwargs.get("simplex_inserted", None) # tried to insert
        self.face             = kwargs.get("face", None) # but this face was not present already
        if self.face and self.simplex_inserted:
            self.message = "You are trying to insert a higher dimensional simplex: {} " \
                           "but one of its faces: {} is not in the Simplicial Complex.".format(
                               self.simplex_inserted, self.face)
    def __str__(self):
        return self.message


class ASC:
    ''' Representation of an Abstract Simplicial Complex in a tree like format
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.root           = Vertex(label = "root", level = -1)
        ## TO DO: vertex_tracker does not seem very useful at the moment. Consider incorporating it into Vertex by making it a dict
        self.vertex_tracker = {}
        self.level_tracker  = collections.defaultdict(set)
        self.max_level      = 0

        ## Init root with some vertices
        for vertex in kwargs["vertices"]:
            vertex.level = 0
            self._add_link(self.root, vertex)
        self.level_tracker[0] = set((self.root,))


    def _add_link(self, parent, vertex):
        ''' Add a connection between parent and child vertices
        '''
        vertex.parent = parent
        parent.connections.add(vertex)
        # track this vertex by some key - not sure if this is the best way
        self.vertex_tracker[(parent.label, vertex.label, vertex.level)] = vertex



    # TO DO: Is order of connections important ? - still a little unsure about how to go about it.
    ## Working in Z_2 for now. But since we are encoding the connections in a tree, a consistent ordering is needed.
    ## Some options to consider:
    ## 1. Make the assumption that labels is sorted by caller in some consistent fashion. (e.g. lexicographically)
    ## 2. Always sort labels and abstract this away from person calling it ?
    ## 3. Provide a pre_sorted option ?
    ## Choosing 2. for now to keep it simple.
    def add_connections(self, labels):
        '''
        :param tuple: labels: vertices having a connection
        '''
        labels = sorted(labels)
        parent = self.vertex_tracker.get(("root", labels[0], 0), None)
        if parent is None:
            raise ASCInsertionError(tuple(labels), tuple((parent.label,)))

        for i in range(0, len(labels) - 1):
            if ("root", labels[i+1], 0) not in self.vertex_tracker:
                raise ASCInsertionError(tuple(labels), tuple((parent.label,)))
            vertex = self.vertex_tracker.get((parent.label, labels[i + 1], i + 1), None)
            if vertex is None:
                ## We always expect a lower dim simplex to be present in this case. E.g. inserting {"Cow", "Horse", "Rabbit"} => {"Cow", "Horse"} is present
                if len(labels) > 2 and (parent.label, labels[i + 1], i) not in self.vertex_tracker:
                    raise ASCInsertionError(tuple(labels), tuple((parent.label, labels[i+1])))
                vertex = Vertex(label = labels[i + 1], parent = parent, level = i + 1)
            self._add_link(parent, vertex)
            parent = vertex

        # keep track of topmost vertex and its levels - might be useful for retreiving n-simplices - NOT SURE
        self.level_tracker[i + 1].add(self.vertex_tracker[("root", labels[0], 0)])
        self.max_level = max(self.max_level, i + 1)


    def _get_simplices(self, n, simplices, parent = None, labels = []):
        ''' recursive method to retrieve labels associated with all n-simplices under a parent
        :param int: n: dimension of the simplex
        :param list: simplices: container for storing any n-simplex found
        :param Vertex: parent: Parent vertex to start the search under
        :param list: labels: ordered list of labels associated with each vertex (top -> bottom)
        '''
        if parent is None:
            parent = self.root

        # add this sorted() iteration, if you want to debug - maybe add a debug flag or something
        # sorted(parent.connections, key = lambda e: e.label):
        for child in parent.connections:
            labels.append(child.label)
            if child.level == n:
                simplices.append(set(labels))
                labels.pop()
            elif child.level > n:
                break
            else:
                self._get_simplices(n, simplices, child, labels)

        if labels:
            labels.pop()


    ## TO DO: extend to allow filter by label/object
    def ret_all_simplices(self, n):
        ''' public method to return all n-simplices
        :param int: n: Return labels associated with all simplices of dimension n in this ASC
        '''
        simplices = []
        seen      = set()
        for i in range(n, self.max_level + 1):
            for parent in self.level_tracker[i]:
                if parent.label in seen:
                    continue
                seen.add(parent.label)
                labels = [parent.label] if parent.label != "root" else []
                self._get_simplices(n, simplices, parent, labels)

        return simplices


    def ret_boundary_matrix(self, dim):
        ''' public method to return the boundary matrix for a particular boundary map
        :param int: dim: corresponding boundary map, e.g. dim = 1 for delta_1; dim = 2 for delta_2, etc.
        '''
        # Get all simplices in C_{dim} and C_{dim-1}
        generators_C_dim         = self.ret_all_simplices(dim)
        generators_C_dim_minus_1 = self.ret_all_simplices(dim - 1) if dim >= 1 \
                                   else [{"0"}]
        if len(generators_C_dim) == 0:
            return np.array([[0]], dtype = np.uint8)
        ## initialize a matrix of zeros - TO DO: Use Z_2/bool for encoding the boundary map/matrix, could not find a built in numpy function to matrix multiply in Z2
        boundary_matrix = np.zeros((len(generators_C_dim_minus_1),
                                    len(generators_C_dim)), dtype = np.uint8)
        ## fill matrix with a 1 if one generator is a boundary of the other
        ## e.g. {Dog} subset {Dog, Horse}; {Dog, Horse} subset {Dog, Horse, Cat}
        for i, r in enumerate(generators_C_dim_minus_1):
            for j, c in enumerate(generators_C_dim):
                if r.issubset(c):
                    boundary_matrix[i][j] = 1

        return (boundary_matrix, generators_C_dim_minus_1, generators_C_dim)







######################################################################################################################################
## Some ideas/notes for future.

## 1.
## Abstract a Simplex as a tree like it is being done in the ASC class and consider the Abstract Simplicial Complex as a Graph
## with each Node as the level 0 node of a Simplex (i.e. a point/vertex) and Edges between them as connections between a point/vertex to another simplex.
## e.g. The 2-simplex A->B->C is say, SimplexX and stored as a tree with level 0 at A;
##      The 1-simplex B->C is say, SimplexY and also stored as a tree with level 0 at B.
## Then the Abstract Simplicial Complex is a Graph like : [SimplexY]---B-@-level0---B-@-level1--->[SimplexX]
## I am thinking only a directed edge from the level0 node to some levelN node is needed. Since the faces are captured in the Simplex itself.
## This allows for queries like : Give me all n-simplices which have the object/label B in it.

## Action to take: Break away tree representation of a Simplex from ASC into a seperate class:
#class Simplex:
#    pass


