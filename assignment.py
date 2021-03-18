## stdlib
import itertools
import sys

## our modules
import abst_simplcl_cmplx

## 3rd party
import numpy as np
import scipy
import scipy.linalg


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


## Took this code as is from the top answer : https://stackoverflow.com/questions/17129290/numpy-2d-and-1d-array-to-latex-bmatrix
def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}_{' + '{} \\times {}'.format(a.shape[0], a.shape[1]) + '}']
    return '\n'.join(rv)


def submission3(data, label):
    ''' Submission for coding assignment 3 - Boundary Matrices
    '''
    ## Using numpy as it prints larger matrices better
    print("Boundary Matrices for {}:\n".format(label))
    for d in [2, 1, 0]:
        ## Get the boundary matrix and print it
        matrix, c_n_1_generators, c_n_generators = data.simplicial_complex.ret_boundary_matrix(d)
        print("$\\delta_{}:".format(d))
        print(bmatrix(matrix))
        print("$")
        print("\n")
        ## store the generators as a tuple(dim > 1) or str(dim==1) for cleaner printing
        c_n_1_generators = [tuple(sorted(e)) if len(e) > 1 else str(next(iter(e))) for e in c_n_1_generators]
        c_n_generators   = [tuple(sorted(e)) if len(e) > 1 else str(next(iter(e))) for e in c_n_generators]
        c_n_vecs = {}
        ## Get all generators in vector form
        for idx, generator in enumerate(c_n_generators):
            c_n_vecs[generator] = np.zeros((len(c_n_generators), 1), dtype = np.bool)
            c_n_vecs[generator][idx] = 1
        ## Boundaries for some linear combinations of our generators
        print("Some select examples of boundaries computed using this matrix :\n")
        for k in range(1, 4):
            for linear_combination in itertools.combinations(c_n_generators, k):
                vector = np.zeros((len(c_n_generators), 1), dtype = np.bool)
                for generator in linear_combination:
                    vector = vector ^ c_n_vecs[generator]

                boundary = (matrix @ vector)%2
                boundary_text = " + ".join(
                    [str(e) for e, mask in zip(c_n_1_generators, boundary) if mask])
                if boundary_text == "":
                    boundary_text = "0"
                linear_combination_text = " + ".join((str(e).replace("'", "") for e in linear_combination))
                print("Boundary of {} is: {}".format(linear_combination_text, boundary_text.replace("'", "")))
                print("\n")
                ## There are too many combinations so just print 1 for each k > 1
                if k > 1:
                    break
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
