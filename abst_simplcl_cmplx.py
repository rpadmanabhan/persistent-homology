## stdlib
import collections
import copy


class Vertex:
    ''' Represents a vertex in an ASC
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.label            = kwargs["label"]
        self.parent           = kwargs.get("parent", None)
        self.level            = kwargs.get("level", None)
        self.connections      = kwargs.get("connections", set()) # child connections
        self.tree_connections = kwargs.get("tree_connections", set()) # connections to edges in another part of the tree

    def __repr__(self):
        '''
        '''
        return "Vertex(label = {}, parent = {}, level = {})".format(self.label, self.parent, self.level)


class ASC:
    ''' Representation of an Abstract Simplicial Complex in a tree like format
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.root           = Vertex(label = "root", level = -1)
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
        # mark a connection to this child vertex from the same labeled vertex at the root level., i.e. level0 to leveln connection
        if parent.label != "root":
            self.vertex_tracker[("root", vertex.label, 0)].tree_connections.add(vertex)
        # track this vertex by some key - not sure if this is the best way
        self.vertex_tracker[(parent.label, vertex.label, vertex.level)] = vertex



    # TO DO: Is order of connections important ?
    def add_connections(self, labels):
        '''
        :param tuple: labels: vertices having a connection
        '''
        parent = self.vertex_tracker[("root", labels[0], 0)]

        for i in range(0, len(labels) - 1):
            vertex = self.vertex_tracker.get((parent.label, labels[i + 1], i + 1), None)
            if vertex is None:
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

        # can take the sorted() out, keeping it for easier debugging
        for child in sorted(parent.connections, key = lambda e: e.label):
            labels.append(child.label)
            if child.level == n:
                simplices.append(list(labels))
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


