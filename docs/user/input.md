# Input File

FUSS reads its settings from `input.ini` in the directory it runs in. The file is parsed by FiNeR.

## Syntax

```ini
[FUSS-Parameters]
res-threshold = 1e-8

[FUSS-Numerics]
time-scheme = RK3
vnn = 2.0               ; inline comment after a value
time-accurate = false
integration-variables = cons
irs-beta = 0.5
```

- A section starts with `[name]`; each option is `key = value`.
- Lines starting with `!`, `;` or `#` are comments, and `;` starts an inline comment after a value.
- FiNeR rejects a comment line that contains `=`: write such notes as an inline comment after a value, or reword them.

---

## Sections Read by FUSS

| Section | Content | Keys |
|---------|---------|------|
| `[FUSS-Parameters]` | New run or restart, stop criteria | [registry](registry.md#fuss-parameters) |
| `[FUSS-Numerics]` | Time scheme, VNN, IRS, integration variables | [registry](registry.md#fuss-numerics) |
| `[FUSS-Multigrid]` | Grid levels and iterations per level | [registry](registry.md#fuss-multigrid) |
| `[FUSS-IO]` | Formats and frequencies of the output | [registry](registry.md#fuss-io) |
| `[FUSS-Probes]` and one section per probe | Probes | [registry](registry.md#fuss-probes), [Running FUSS](using.md#probes) |

The [Parameter Registry](registry.md) lists every key with its default and allowed values. It is generated from the source by `bin/DocGen` (`src/app/docgen.f90`): do not edit it by hand. Keys that come in numbered series appear there once as examples: `probe<n>` with its section (`probe-placeholder`), and `level<m>-iter` for two levels.

One key is not in the registry: `levels` in `[FUSS-Multigrid]`, the number of grid levels (default 1, no multigrid). FUSS reads it before building the registry, to know how many `level<m>-iter` keys to expect. Multigrid is described in [Time Integration](../theory/time-integration.md#multigrid-grid-sequencing).

The other sections of the test-case files are not read by FUSS; they are the input of the pre-processors that write `INPUT/`:

| Section | Pre-processor | Described in |
|---------|---------------|--------------|
| `[GRIB-meshgen]`, `[GRIB-Block<n>]` | GRIB (mesh) | [Running FUSS](using.md#mesh) |
| `[GPB-Phase1]` | ATLAS GPB (`phase.txt`, `properties.dat`) | [Initial Conditions](initial-conditions.md#materials) |
| `[ICB-Block<n>]` | ATLAS ICB (`ic.tec`) | [Initial Conditions](initial-conditions.md#initial-field-inputic) |
| `[STB-Block<n>]` | ATLAS STB (`st.tec`) | [Initial Conditions](initial-conditions.md#heat-source-inputst) |
| `[ATLAS-Parameters]`, `[BCB-Block<n>]`, named BC sections | ATLAS BCB (`bc.txt`) | [Boundary Conditions](boundary-conditions.md#atlas-input) |

---

## Defaults and Validation

- A key that is not set takes its default, except `time-accurate`, which is required: without it the run stops at set-up.
- FUSS looks up only the keys it knows: a misspelled key or section is ignored without a message, and the default is used. Check the spelling against the [registry](registry.md).
- Values are checked at set-up against the allowed values and ranges of the registry. A value out of range stops the run with `[ERROR] <key> must be ...` under `Input file FAIL`. Text values are case-sensitive: `time-scheme = rk3` is rejected.
- IRS is switched on by `irs-beta > 0`.
- With more than one multigrid level `iter-threshold` has no effect: each level stops at its `level<m>-iter`. A run that sets it anyway gets a `WARNING` line under `Input file OK`.

## Changes During a Run

FUSS reads `input.ini` again every `ini-diter` iterations (default 10000). Keys used during the iterations take effect from then on, for example `vnn`, `res-threshold`, `time-threshold`, `iter-threshold`, `level<m>-iter` and the output frequencies. IRS cannot be switched on at runtime: `irs-beta > 0` in a run started without it is ignored with a warning, while `irs-beta = 0` switches it off (see [Time Integration](../theory/time-integration.md#implicit-residual-smoothing-irs)). The values read at runtime are not checked against the registry.
