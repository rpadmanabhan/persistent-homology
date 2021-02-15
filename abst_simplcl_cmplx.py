## stdlib
import collections
import copy


class Vertex:
    ''' Represents a vertex in an ASC
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.label       = kwargs["label"]
        self.parent      = kwargs.get("parent", None)
        self.level       = kwargs.get("level", None)
        self.connections = kwargs.get("connections", set())

    def __repr__(self):
        '''
        '''
        return "Vertex(label = {}, parent = {}, level = {})".format(self.label, self.parent, self.level)


class ASC:
    ''' Abstract Simplicial Complex
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.root           = Vertex(label = "root", level = -1)
        self.vertex_tracker = {}

        ## Don't use these for anything at the moment
        self.level_tracker  = collections.defaultdict(set)
        self.max_level      = 0

        ## Init root with some vertices
        for vertex in kwargs["vertices"]:
            vertex.level = 0
            self._add_link(self.root, vertex)


    def _add_link(self, parent, vertex):
        ''' Add a connection between parent and child vertices
        '''
        vertex.parent = parent
        parent.connections.add(vertex)
        self.vertex_tracker[(parent.label, vertex.label, vertex.level)] = vertex


    def _ret_simplices(self, n, simplices, parent = None, labels = []):
        ''' internal method to traverse tree structure and return n-simplexes
        '''
        if parent is None:
            parent = self.root

        for child in parent.connections:
            labels.append(child.label)
            if len(labels[child.level - n:]) == n + 1:
                simplices.append(labels[child.level - n:])
            if child.level >= n:
                self._ret_simplices(n, simplices, child, labels[n - child.level + 1:])
            else:
                self._ret_simplices(n, simplices, child, labels)

            labels.pop()


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


    def ret_all_simplices(self, n):
        ''' public method to return all n-simplices
        :param int: n: Return labels associated with all simplices of dimension n in this ASC
        '''
        simplices = []
        self._ret_simplices(n, simplices)
        return simplices

