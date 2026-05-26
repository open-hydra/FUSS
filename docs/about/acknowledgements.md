# Acknowledgements

FUSS is built upon several open-source projects.

## FiNeR

**Fortran INI ParseR and generator**

- **Repository:** [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR)
- **License:** GPL v3.0

FiNeR is a pure Fortran 2003+ OOP library for reading and writing INI configuration files. FUSS uses FiNeR to parse the `input.ini` parameter file.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read and write files in different formats. FUSS uses ORION for solution output in Tecplot and VTK formats and for reading the structured-grid mesh.

## Documentation Tools

### MkDocs

**Static Site Generator**

- **Website:** [mkdocs.org](https://www.mkdocs.org)
- **License:** BSD-2-Clause
- **Contribution:** Documentation site generation

MkDocs transforms FUSS's documentation into a beautiful, searchable website.

### Material for MkDocs

**Modern Documentation Theme**

- **Website:** [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
- **License:** MIT
- **Contribution:** Professional documentation theme with advanced features

Material for MkDocs provides the sleek, modern interface for FUSS's documentation.

## Institutional Support

While FUSS is an independent open-source project, it has benefited from:

- Academic research environments
- Access to the NAFEMS thermal benchmark suite and other reference datasets for validation

## License Compliance

FUSS respects all licenses of dependencies and foundations:

- **GPL v3.0** — ORION, FiNeR

See the [License](license.md) page for FUSS's full license text.

---