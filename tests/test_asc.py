## stdlib
import os
import sys
import unittest

## our modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import abst_simplcl_cmplx



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

        self.assertEqual(self.simplicial_complex.vertex_tracker[("Rabbit", "Horse", 2)].level, 2)
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Cow", "Horse", 1)].level, 1)
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Rabbit", "Horse", 1)].level, 1)


        self.assertEqual(self.simplicial_complex.ret_all_simplices(2), [["Cow", "Rabbit", "Horse"]])
        self.assertEqual(sorted(self.simplicial_complex.ret_all_simplices(1)), sorted([["Cow", "Rabbit"], ["Cow", "Horse"], ["Rabbit", "Horse"]]))
        self.assertEqual(sorted(self.simplicial_complex.ret_all_simplices(0)), sorted([["Cow"], ["Rabbit"], ["Horse"], ["Dog"]]))


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
        '''
        '''
        vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in \
             ["Cow", "Rabbit", "Horse", "Dog", "Fish",
              "Dolphin", "Oyster", "Broccoli", "Fern", "Onion", "Apple"]))
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = vertices)
        for connections in [
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
                ("Fern", "Onion", "Apple")]:
            simplicial_complex.add_connections(connections)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(2)), 9)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(1)), 15)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(0)), len(vertices))


    def test_exampleB(self):
        '''
        '''
        vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in \
             ["Cow", "Rabbit", "Horse", "Dog", "Fish",
              "Dolphin", "Oyster", "Broccoli", "Fern", "Onion", "Apple"]))
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = vertices)
        for connections in [
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
                ("Oyster", "Broccoli", "Apple"), ("Oyster", "Onion", "Apple")]:
            simplicial_complex.add_connections(connections)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(2)), 8)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(1)), 27)
        self.assertEqual(len(simplicial_complex.ret_all_simplices(0)), len(vertices))



if __name__ == '__main__':
    unittest.main()
