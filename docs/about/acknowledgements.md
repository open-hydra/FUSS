# Acknowledgements

FUSS is built upon several open-source projects.

## FiNeR

**Fortran INI ParseR and generator**

- **Repository:** [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR)
- **License:** multi-licensed: GPL v3.0, BSD 2-Clause, BSD 3-Clause or MIT

FiNeR is a pure Fortran 2003+ OOP library for reading and writing INI configuration files. FUSS uses FiNeR to parse the `input.ini` parameter file. FiNeR depends on PENF, StringiFor, FACE, FLAP and BeFoR64 (same author, same licenses), which `install.sh` clones.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read and write files in different formats. FUSS uses ORION to read the initial condition (grid and fields), the property tables and restart files, and to write the solution, wall and diagnostic files in Tecplot and VTK formats.

## TecIO

**Tecplot binary I/O library** (optional)

- **Provider:** Tecplot, Inc.
- **License:** Tecplot's TecIO license agreement (`lib/TecIO/teciosrc/tecio_license_agreement.txt` in ORION)

TecIO ships with ORION and is built only with `--use-tecio`; FUSS then writes Tecplot binary (`.szplt`) files.

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

- **GPL v3.0** — ORION
- **GPL v3.0 / BSD / MIT** (multi-licensed) — FiNeR and its dependencies
- **TecIO license agreement** — TecIO (optional)

See the [License](license.md) page for FUSS's full license text.

---