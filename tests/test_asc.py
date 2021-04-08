## stdlib
import os
import pprint
import sys
import unittest

## our modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import abst_simplcl_cmplx
import assignment

## 3rd party
import numpy as np

# Team : Bill Lee, Raghavendra Padmanabhan and Francisco Vargas
#--------------------------------------------------------------------------------------

class TestASCInit(unittest.TestCase):

    def setUp(self):
        '''
        '''
        self.vertices = set((abst_simplcl_cmplx.Vertex(label = "Cow"),
                            abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                            abst_simplcl_cmplx.Vertex(label = "Horse"),
                            abst_simplcl_cmplx.Vertex(label = "Dog")))

    def test_init_ASC(self):
        ''' Test if vertices are correctly initialized
        '''
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)
        self.assertEqual(
            simplicial_complex.root.connections,
            self.vertices
        )

    def tearDown(self):
        '''
        '''
        pass


class TestASCInsertConnection(unittest.TestCase):

    def setUp(self):
        '''
        '''
        self.vertices = set((abst_simplcl_cmplx.Vertex(label = "Cow"),
                             abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                             abst_simplcl_cmplx.Vertex(label = "Horse"),
                             abst_simplcl_cmplx.Vertex(label = "Dog")))
        self.simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)


    def test_insert_connections(self):
        ''' Insert some connections and test if correct structure is maintained
        '''
        self.simplicial_complex.add_connections(("Cow", "Rabbit"))
        self.simplicial_complex.add_connections(("Cow", "Horse"))
        self.simplicial_complex.add_connections(("Rabbit", "Horse"))
        self.simplicial_complex.add_connections(("Cow", "Rabbit", "Horse"))
        ## Note: order is different here, because add_connections always sorts lexicographically.
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Horse", "Rabbit", 2)].level, 2)
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Cow", "Horse", 1)].level, 1)
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Horse", "Rabbit", 1)].level, 1)

        self.assertEqual(self.simplicial_complex.ret_all_simplices(2),
                         [{"Cow", "Rabbit", "Horse"}])
        self.assertEqual(sorted(self.simplicial_complex.ret_all_simplices(1), key = sorted),
                         sorted([{"Cow", "Rabbit"}, {"Cow", "Horse"}, {"Rabbit", "Horse"}], key = sorted))
        self.assertEqual(sorted(self.simplicial_complex.ret_all_simplices(0), key = sorted),
                         sorted([{"Cow"}, {"Rabbit"}, {"Horse"}, {"Dog"}], key = sorted))


    def test_insert_connections_invalid(self):
        ''' Make sure we throw an exception when an invalid insertion takes place.
        '''
        self.assertRaises(abst_simplcl_cmplx.ASCInsertionError,
                          self.simplicial_complex.add_connections,
                          ("Zebra", "Cow"))
        self.simplicial_complex.add_connections(("Cow", "Rabbit"))
        self.assertRaises(abst_simplcl_cmplx.ASCInsertionError,
                          self.simplicial_complex.add_connections,
                          ("Cow", "Rabbit", "Horse"))
        self.simplicial_complex.add_connections(("Cow", "Horse"))
        self.assertRaises(abst_simplcl_cmplx.ASCInsertionError,
                          self.simplicial_complex.add_connections,
                          ("Cow", "Rabbit", "Horse"))



    def tearDown(self):
        '''
        '''
        pass


class TestExamples(unittest.TestCase):
    def setUp(self):
        '''
        '''
        pass

    def test_exampleA(self):
        ''' Part A data from Coding assignment 1
        '''
        vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in \
             assignment.partA.vertex_labels))
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = vertices)
        for connections in assignment.partA.vertex_connections:
            simplicial_complex.add_connections(connections)

        ## Check we have exactly the same simplices we inserted
        self.assertEqual(sorted(simplicial_complex.ret_all_simplices(2), key = sorted),
                         sorted([set(conn) for conn in assignment.partA.vertex_connections \
                                 if len(conn) == 3], key = sorted))
        self.assertEqual(sorted(simplicial_complex.ret_all_simplices(1), key = sorted),
                         sorted([set(conn) for conn in assignment.partA.vertex_connections \
                                 if len(conn) == 2], key = sorted))
        self.assertEqual(sorted(simplicial_complex.ret_all_simplices(0), key = sorted),
                         sorted([set((e,)) for e in assignment.partA.vertex_labels], key = sorted))
        ## Double check that we have the expected number of them
        self.assertEqual(len(simplicial_complex.ret_all_simplices(2)), 9)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(1)), 15)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(0)), len(vertices))

    def test_exampleB(self):
        ''' Part B data from Coding assignment 1
        '''
        vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in \
             assignment.partB.vertex_labels))
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = vertices)
        for connections in assignment.partB.vertex_connections:
            simplicial_complex.add_connections(connections)
        self.assertEqual(sorted(simplicial_complex.ret_all_simplices(2), key = sorted),
                         sorted([set(conn) for conn in assignment.partB.vertex_connections \
                                 if len(conn) == 3], key = sorted))
        self.assertEqual(sorted(simplicial_complex.ret_all_simplices(1), key = sorted),
                         sorted([set(conn) for conn in assignment.partB.vertex_connections \
                                 if len(conn) == 2], key = sorted))
        self.assertEqual(sorted(simplicial_complex.ret_all_simplices(0), key = sorted),
                         sorted([set((e,)) for e in assignment.partB.vertex_labels], key = sorted))
        ## Double check that we have the expected number of them
        self.assertEqual(len(simplicial_complex.ret_all_simplices(2)), 8)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(1)), 27)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(0)), len(vertices))


class RetBoundaryMatrix(unittest.TestCase):
    def setUp(self):
        '''
        '''
        pass

    def test_2simplex(self):
        ''' Compute boundary matrices for different maps of a 2-simplex, i.e. a triangle (filled)
        '''
        vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) \
             for label in ["V1", "V2", "V3"]))
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = vertices)
        for connections in [("V1", "V2"), ("V1", "V3"), ("V2", "V3"),
                            ("V1", "V2", "V3")]:
            simplicial_complex.add_connections(connections)
        self.assertEqual(
            (simplicial_complex.ret_boundary_matrix(0)[0] ==  \
             np.array([[0, 0, 0]], dtype = int)).all(), True)
        # Note: Order can be different depending on how data was loaded to tree
        # sorted is only to rearrange columns of matrix for comparison for equality - no change in interpretation of the matrix.
        self.assertEqual(
            (np.sort(simplicial_complex.ret_boundary_matrix(1)[0]) == \
             np.sort(np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]],
                              dtype = int))).all(), True)
        self.assertEqual(
            (simplicial_complex.ret_boundary_matrix(2)[0] == \
             np.array([[1], [1], [1]], dtype  = int)).all(), True)
        self.assertEqual(
            (simplicial_complex.ret_boundary_matrix(3)[0] == \
             np.array([[0]], dtype = int)).all(), True)


    def test_boundary_of_boundary(self):
        ''' Test to make sure the boundary of a boundary is 0
        '''
        vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in \
             assignment.partA.vertex_labels))
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = vertices)
        for connections in assignment.partA.vertex_connections:
            simplicial_complex.add_connections(connections)

        boundary_matrix_C_2, c_1_generators, c_2_generators = simplicial_complex.ret_boundary_matrix(2)
        boundary_matrix_C_1, c_0_generators, _ = simplicial_complex.ret_boundary_matrix(1)
        for _ in range(100):
            rand_vec = np.random.randint(2, size = (len(c_2_generators), 1), dtype = np.uint8)
            boundary1 = (boundary_matrix_C_2 @ rand_vec) % 2
            boundary2 = (boundary_matrix_C_1 @ boundary1) % 2
            self.assertEqual(np.all(boundary2 == 0), True)

        boundary_matrix_C_0, _, _ = simplicial_complex.ret_boundary_matrix(0)
        for _ in range(100):
            rand_vec = np.random.randint(2, size = (len(c_1_generators), 1), dtype = np.uint8)
            boundary1 = (boundary_matrix_C_1 @ rand_vec) % 2
            boundary2 = (boundary_matrix_C_0 @ boundary1) % 2
            self.assertEqual(np.all(boundary2 == 0), True)

if __name__ == '__main__':
    unittest.main()
