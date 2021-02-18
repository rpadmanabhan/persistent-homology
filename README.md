# abstract-simplicial-complex

## Representation of an abstract simplicial complex as a tree datastructure in Python3.

#### Usage - Assuming python3.5+, should work for any python3
```python
import abst_simplcl_cmplx

## Initialize an abstract simplicial complex
simplicial_complex = abst_simplcl_cmplx.ASC(
          vertices = [abst_simplcl_cmplx.Vertex(label = "Cow"), abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                      abst_simplcl_cmplx.Vertex(label = "Rabbit")])
## Add some connections
simplicial_complex.add_connections(("Cow", "Rabbit"))
simplicial_complex.add_connections(("Cow", "Horse"))
simplicial_complex.add_connections(("Rabbit", "Horse"))
simplicial_complex.add_connections(("Cow", "Rabbit", "Horse"))

## Return back the 2-simplex just inserted above
simplicial_complex.ret_all_simplices(2)
```

#### Run tests
```
python3 tests/test_asc.py -v
```

##### Note: abstract simplicial complexes A and B from the coding assignment are represented in class TestExamples under tests/test_asc.py

