## stdlib
import pprint

## our modules
import abst_simplcl_cmplx

## 3rd part
import numpy as np


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
    ''' Submission for coding assignment 2
    '''
    print(f"Boundary Matrix for {label}:\n")
    print("del2:")
    pprint.pprint(np.matrix(data.simplicial_complex.ret_boundary_matrix(2)))
    print("del1:")
    pprint.pprint(np.matrix(data.simplicial_complex.ret_boundary_matrix(1)))
    print("del0:")
    pprint.pprint(np.matrix(data.simplicial_complex.ret_boundary_matrix(0)))
    print("\n")


if __name__ == '__main__':
    submission2(partA, "Part A")
    submission2(partB, "Part B")
