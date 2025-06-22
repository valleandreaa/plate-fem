[MESHER]
- able to generate the mesh and coverting it as necessary.
[FEM-Engine]
- applying physical parameters and forcess through a terminal. which also opens the mesh assigning references in order to help out in the forces assignment

[SOLVER]
- being able to apply different solvers 

[OUTPUT-FILE]
- being able to output the files

[OVERALL]
Make a script that through a terminal, is able to:
- assign the nodes
- then the segmentes or arc
- then physical properties/thickness if necessary
- then boundaries and forces
- automatic meshing (choose paraemters and show finess)
- select solvers
- print out the files to visualize results.



PlateFEM provides simple finite element utilities for plate analysis. A small interactive CLI allows creation of nodes, segments and arcs, definition of material parameters and solution of the problem with a linear solver.

The library lives in the `platefem` package. Documentation can be generated with Sphinx:

```bash
sphinx-build -b html docs/source docs/build/html
```

The CLI entry point can be executed with:

```bash
python -m platefem.cli
```