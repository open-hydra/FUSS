# Versioning

FUSS requires ORION by a range of releases and FiNeR at one release. Both are checked at configure
time, so a library FUSS was not tested with stops the configure before anything is built.

## How FUSS requires its libraries

`CMakeLists.txt` names the releases FUSS is tested with:

| Library | Requirement | How |
|---|---|---|
| ORION | `FUSS_ORION_VERSION 2.1`: 2.1.0 and later, up to 3.0.0 excluded | `semver_dependency` of `cmake/SemVerDeps.cmake` |
| FiNeR | `FUSS_FINER_TAG v2.3.0` | the submodule `lib/FiNeR`, checked against its `VERSION` file |

`cmake/SemVerDeps.cmake` is FLINT's file, copied unmodified: update it from FLINT, never here.
A requirement accepts the releases from its own up to the next change of its leftmost nonzero
number (`2.1` accepts 2.1.0 to 3.0.0 excluded, `0.3` accepts 0.3.x only). Only `vX.Y.Z` tags are
releases; a tag such as `v2.2.0-rc1` is never selected.

FUSS requires the oldest release that has everything it uses:

- **ORION `2.1`**: the Tecplot ASCII reader that FUSS uses for `ic.tec`, `st.tec` and the restart
  solutions reads data lines of any length, and stops at an error in the coordinates of a BLOCK
  zone.
- **FiNeR**: FUSS needs v2.1.0 or later, since an older FiNeR ends the reading of a section of
  `input.ini` at a comment line that holds `=`. It pins v2.3.0.

FiNeR's own dependencies (PENF, StringiFor, FACE, FLAP, BeFoR64) are fetched as FiNeR's author
does (`fobis fetch`: their default branches, `FUSS_FINER_DEPS=FREE`). If an upstream change breaks
FUSS, `FUSS_FINER_DEPS=BACKUP` checks out the releases of `FUSS_FINER_BACKUP_PINS`.

## The checks

| Library | What stops the configure |
|---|---|
| `lib/FiNeR` (submodule) | checked out at another commit than the gitlink of this FUSS commit (after a `git checkout` without `git submodule update`), or its `VERSION` differs from `FUSS_FINER_TAG` |
| FiNeR clone (`--include-finer`) | its `VERSION` differs from `FUSS_FINER_TAG`; with `BACKUP`, a dependency not at its release |
| ORION clone (`--include-orion`) | its newest `vX.Y.Z` tag is not in the required range |

To build against a work-in-progress checkout on purpose (an ORION or FiNeR branch under test),
configure with `-DFUSS_PIN_CHECK=WARNING`: each check then warns instead.

## Which release a build gets

The release of ORION is, in order:

1. the clone given with `-DFETCHCONTENT_SOURCE_DIR_ORION=<folder>` (`install.sh --include-orion`),
   as it is; the configure stops if it is not a release that FUSS accepts;
2. the release in `deps.lock`, if it satisfies the requirement;
3. the newest release on GitHub that satisfies it.

The configure checks the release out in `lib/ORION` and writes it to `deps.lock`; both are git
ignored. `-DSEMVER_UPDATE=ALL` (`install.sh build --update-deps`) selects the newest one again.

FiNeR is the submodule `lib/FiNeR`, or the clone given with `-DFETCHCONTENT_SOURCE_DIR_FINER`
(`install.sh --include-finer`). The configure clones FiNeR's dependencies into
`lib/FiNeR/src/third_party/`, ignored by FiNeR's repository.

## The requirements in hydra

Hydra builds FUSS with its own ORION and FiNeR. When hydra selects ORION with SemVerDeps, FUSS's
requirement is checked against that release, so a new major ORION in FUSS's requirement, or another
FiNeR release, must move in hydra too.

## Bumping a library

- ORION: raise `FUSS_ORION_VERSION` in `CMakeLists.txt`, configure with `--update-deps`, run the
  suite.
- FiNeR: check the new release out in `lib/FiNeR`, `git add lib/FiNeR`, set `FUSS_FINER_TAG` to it
  in the same commit, run the suite.
- FiNeR's dependencies under `BACKUP`: edit `FUSS_FINER_BACKUP_PINS` in `CMakeLists.txt`.
