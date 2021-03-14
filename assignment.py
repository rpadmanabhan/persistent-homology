## stdlib
import itertools
import pprint
import sys

## our modules
import abst_simplcl_cmplx

## 3rd party
import numpy as np



# Team : Bill Lee, Raghavendra Padmanabhan and Francisco Vargas
#--------------------------------------------------------------------------------------


class DataCoding1:
    ''' Simplicial Complex for the data in Coding assignment 1
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.vertex_labels      = kwargs["vertex_labels"]
        self.vertex_connections = kwargs["vertex_connections"]
        ## Build the simplicial complex
        self.vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in self.vertex_labels))
        self.simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)
        for connections in self.vertex_connections:
            self.simplicial_complex.add_connections(connections)

## Part A
partA = DataCoding1(
    vertex_labels = set(
        {"Cow", "Rabbit", "Horse", "Dog", "Fish",
         "Dolphin", "Oyster", "Broccoli", "Fern", "Onion", "Apple"}),
    vertex_connections = [
        ("Cow", "Rabbit"), ("Cow", "Horse"), ("Cow", "Dog"),
        ("Rabbit", "Horse"), ("Rabbit", "Dog"), ("Horse", "Dog"),
        ("Fish", "Dolphin"), ("Fish", "Oyster"), ("Dolphin", "Oyster"),
        ("Broccoli", "Fern"), ("Broccoli", "Onion"), ("Broccoli", "Apple"),
        ("Fern", "Onion"), ("Fern", "Apple"),
        ("Onion", "Apple"),
        ("Cow", "Rabbit", "Horse"), ("Cow", "Rabbit", "Dog"),
        ("Cow", "Horse", "Dog"), ("Rabbit", "Horse", "Dog"),
        ("Fish", "Dolphin", "Oyster"),
        ("Broccoli", "Fern", "Onion"),
        ("Broccoli", "Fern", "Apple"),
        ("Broccoli", "Onion", "Apple"),
        ("Fern", "Onion", "Apple")])

## Part B
partB = DataCoding1(
    vertex_labels = set(
        {"Cow", "Rabbit", "Horse", "Dog", "Fish",
         "Dolphin", "Oyster", "Broccoli", "Fern", "Onion", "Apple"}),
    vertex_connections = [
        ("Cow", "Rabbit"), ("Cow", "Fish"), ("Cow", "Oyster"),
        ("Cow", "Broccoli"), ("Cow", "Onion"), ("Cow", "Apple"),
        ("Rabbit", "Fish"), ("Rabbit", "Oyster"), ("Rabbit", "Broccoli"),
        ("Rabbit", "Onion"), ("Rabbit", "Apple"), ("Fish", "Oyster"),
        ("Fish", "Broccoli"), ("Fish", "Onion"), ("Fish", "Apple"),
        ("Oyster", "Broccoli"), ("Oyster", "Onion"), ("Oyster", "Apple"),
        ("Broccoli", "Onion"), ("Broccoli", "Apple"), ("Onion", "Apple"),
        ("Horse", "Dog"), ("Horse", "Dolphin"), ("Horse", "Fern"),
        ("Dog", "Dolphin"), ("Dog", "Fern"), ("Dolphin", "Fern"),
        ("Cow", "Broccoli", "Apple"), ("Cow", "Onion", "Apple"),
        ("Rabbit", "Broccoli", "Apple"), ("Rabbit", "Onion", "Apple"),
        ("Fish", "Broccoli", "Apple"), ("Fish", "Onion", "Apple"),
        ("Oyster", "Broccoli", "Apple"), ("Oyster", "Onion", "Apple")])



def submission2(data, label):
    ''' Submission for coding assignment 2 - Boundaries of the p-chains.
    '''
    ## Boundaries in C_2 are the boundaries of all 2-simplices, i.e. all edges or 1-simplices
    print("Boundaries in $C_2$ are: \n")
    for simplex in data.simplicial_complex.ret_all_simplices(2):
        out = []
        for b in itertools.combinations(list(simplex), 2):
            out.append(str(tuple(b)))
        out = " + ".join(out)
        print("Boundary of the 2-simplex: {} is: {}\n".format(tuple(simplex), out))
    ## Boundaries in C_1 are the boundaries of all 1-simplices, i.e. all vertices or  0-simplices
    print("\nBoundaries in $C_1$ are: \n")
    for simplex in data.simplicial_complex.ret_all_simplices(1):
        v1, v2 = list(simplex)
        print("Boundary of the 1-simplex: {} is: {} + {}\n".format(tuple(simplex), v1, v2))
    ## Boundaries of C_0 are the boundaries of all 0-simplices, i.e. just 0
    print("\nBoundaries in $C_0$ are: \n")
    print("\nBoundary of a vertex is 0, so all boundaries in $C_0$ are just 0.")

def submission3(data, label):
    ''' Submission for coding assignment 3 - Boundary Matrices
    '''
    ## Using numpy as it prints larger matrices better
    print("Boundary Matrix for {}:\n".format(label))
    print("del2:")
    pprint.pprint(np.matrix(data.simplicial_complex.ret_boundary_matrix(2)))
    print("del1:")
    pprint.pprint(np.matrix(data.simplicial_complex.ret_boundary_matrix(1)))
    print("del0:")
    pprint.pprint(np.matrix(data.simplicial_complex.ret_boundary_matrix(0)))
    print("\n")



if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Specify param coding2 or coding3 for specific output")
        sys.exit(0)
    if sys.argv[1].lower() == "coding2":
        print("PartA\n")
        submission2(partA, "Part A")
        print("\n-------------------------------------------------------------------\n")
        print("PartB\n")
        submission2(partB, "Part B")
        print("\n-------------------------------------------------------------------\n")
    elif sys.argv[1].lower() == "coding3":
        submission3(partA, "Part A")
        submission3(partB, "Part B")
