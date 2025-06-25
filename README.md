# PlateFEM

PlateFEM is a small collection of utilities for 2D plate finite element analysis.
It provides a minimal mesher, a FEM engine and different solver backends that can
be invoked from a command line interface.  The project originated from
university assignments and is organised as a standard Python package.

## Features

- **Mesher** – generates triangular or quadrilateral meshes using Gmsh and
  reports the element skewness.
- **FEM engine** – applies material properties, loads and boundary conditions and
  computes element stresses.
- **Solvers** – linear and non‑linear solvers are provided in the `platefem.solvers`
  module.
- **CLI** – create models from YAML configuration files and write VTK results for
  visualisation.

Example configuration files can be found in `sample_config.yaml` and simple usage
scripts are located in the `examples` folder.

## Sample output

Running the cli for a plate and displaying the VTK results we get:

![Sample Von Mises](examples\imgs\VonMises_8000.png)

## Documentation

The API is documented with Sphinx.  HTML pages can be generated with:

```bash
sphinx-build -b html docs/source docs/build/html
```

## License

This project is distributed for educational purposes and comes without any
warranty.