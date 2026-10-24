#!/usr/bin/env python3
"""
Split the single block of a FUSS case into nbi x nbj blocks (plan 11 stage 2).

    usage: split_blocks.py [--set SECTION:KEY=VALUE]... <case_dir> <nbi> <nbj> <out_dir>

--set replaces the value of KEY in the kept section [SECTION] (the key must be
there): a twin with a fixed cadence, say, where the source stops on its residual
(plan 11 stage 3: a 1- vs N-rank comparison cannot stop on res-threshold).

WHY
---
Every moving-mesh case of the suite is single-block, so running one under 4 MPI
ranks would test one working rank and three idle ones. This tool writes a
multi-block twin of such a case: the same mesh, the same initial field and the
same boundary conditions, cut into blocks joined by 101 connections, so the 1-
vs N-rank comparison of plan 11 has blocks to distribute.

WHAT IT WRITES
--------------
The block is cut at node indices into nbi blocks along i and nbj along j, equal
shares (refused when the cells do not divide evenly); K is not cut. Blocks are
numbered row-major, b = 1 + bi + nbi*bj (bi along i, bj along j).

  INPUT/ic.tec     one zone per block. The VARIABLES line is copied, each zone
                   header is the source header with only the title (B<b>-SP)
                   and I/J/K changed. Nodal x, y, z, then the cell-centred
                   variables, one value per line. The values are the source's
                   own STRINGS, re-arranged, never re-formatted floats, so every
                   value is bit-identical to the unsplit case's.
  MESH/mesh.tec    the same way, titles Block<b>, if the source has one (FUSS
                   does not read it; it is kept consistent with ic.tec). The
                   titles follow test/moving-mesh/conjugate_receding, the
                   suite's multi-block case: B<b>-SP in ic.tec, Block<b> in
                   mesh.tec.
  INPUT/bc.txt     block-major, faces 1..6 in order, each (block, face)
                   contiguous: Mod_BC_Fluxes and IO_Wall walk the records by the
                   cumulative counts n_bf(b,f), so the grouping is load-bearing
                   (Check_BC only counts records). Each block carries
                   2*nj*nk + 2*ni*nk + 2*ni*nj records (cells). Within a face
                   the records follow the order of the SOURCE face (its record
                   sequence restricted to the block's cells), whatever loop the
                   source was written with.
                     outer face: the source record's type and second line,
                                 copied verbatim, re-indexed to the block's cells
                                 (a 307's second line names the coupled solver's
                                 cell, not a FUSS one, so it is kept as it is);
                     cut face:   type 101, second line
                                 bs is js ks fs 1 0 0 1
                                 (the neighbour block's adjacent cell and face,
                                 identity orientation; the form of conn() in
                                 test/numerics-features/conjugate_two_layer/
                                 generate.py, read at IO_BC.f90 Read_BCfile).
  INPUT/phase.txt, INPUT/properties.dat   copied, and the per-case tables of the
                   pyrolysis, Darcy and ablation cases (plan 11 stage 3):
                   pyrolysis.dat, virgin.dat, char.dat, gas.dat, env_*.dat,
                   bprime_*.dat -- one file for every block, nothing per cell --
                   and the TEST ONLY [FUSS-Ablation] coupled-env-file when it is
                   one table for every 307 face (plan 11 stage 4).
  FUSS.sh          copied, MASTERDIR set for out_dir's depth below the
                   repository (../../.. for test/moving-mesh/<case>).
  input.ini        only the [FUSS-*] sections, verbatim, behind a one-line
                   comment naming the source. The ATLAS sections ([GRIB-*],
                   [GPB-*], [ICB-*], [BCB-*]) and the patch sections they name
                   describe the unsplit case and are dropped, with the source's
                   preamble. [FUSS-Probes] is kept with the section of every
                   probe it names, an index-position "1 i j k" re-indexed to
                   the block that holds cell (i, j) (a position by coordinates
                   is kept as it is). No comment line may hold '=' (rule of every
                   hydra INI): a kept one that does is refused, not dropped.

Files the tool does not write in out_dir (verify.py, README.md) are left alone.

WHAT IT REFUSES
---------------
A source with more than one zone; records that already carry connectivity
(101, 102, 103, 104: their targets would need re-mapping); a bc.txt that
does not cover every boundary cell exactly once, face by face in order 1..6;
INPUT files it does not know how to split (st.tec, multigrid bcN.txt, ...);
a coupled-env-file of "// face b i j k f" blocks (they name the source's bc
records); multigrid sections (coarse-level files); time or probe files at the
case root.

Python 3.6, standard library only.
"""

import os
import re
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # test/tools -> repository

# Record layout of bc.txt, as Check_BC (IO_BC.f90) reads it: these types carry
# one value line after "b i j k f type"; 102/104 carry a count line and donor
# lines; every other type (0 null, 300 symmetry) carries none.
ONE_VALUE_LINE = (101, 103, 301, 302, 303, 304, 305, 306, 307)
# (307 is not here: its value line names the coupled solver's cell, not a FUSS one)
CONNECTIVITY = {101: "block connection", 102: "chimera", 103: "coupled wall",
                104: "coupled chimera"}

KNOWN_INPUT = ("ic.tec", "bc.txt", "phase.txt", "properties.dat")
# copied verbatim: tables of the whole case, not of a block or a cell (plan 11 stage 3)
COPY_INPUT = ("phase.txt", "properties.dat", "pyrolysis.dat", "virgin.dat", "char.dat", "gas.dat")
COPY_INPUT_PATTERNS = (r"^env_[\w.-]*\.dat$", r"^bprime_[\w.-]*\.dat$")


COPY_EXTRA = set()                                       # the coupled-env-file (main)


def copied_input(name):
    return name in COPY_INPUT or name in COPY_EXTRA or any(re.match(p, name) for p in COPY_INPUT_PATTERNS)


def coupled_env_file(text, src):
    """The TEST ONLY [FUSS-Ablation] coupled-env-file of input.ini, or None. One table
    for every 307 face is copied; a file of "// face b i j k f" blocks names the
    source's bc records and is refused."""
    sec, name = None, None
    for line in text.split("\n"):
        m = re.match(r"^\s*\[([^\]]+)\]", line)
        if m:
            sec = m.group(1).strip()
            continue
        m = re.match(r"^\s*coupled-env-file\s*=\s*([^;!#\s]+)", line)
        if sec == "FUSS-Ablation" and m:
            name = m.group(1)
    if name is None or name.lower() == "none":
        return None
    path = os.path.join(src, "INPUT", name)
    if not os.path.isfile(path):
        die("coupled-env-file = %s: no INPUT/%s" % (name, name))
    with open(path) as fh:
        if any(re.match(r"^\s*//\s*face\b", l) for l in fh):
            die("INPUT/%s holds '// face' blocks: they name the source's bc records "
                "(not re-mapped)" % name)
    return name
COMMENT_CHARS = "!;#"                                   # FiNeR's COMMENTS


def die(msg):
    sys.stderr.write("split_blocks.py: ERROR: %s\n" % msg)
    sys.exit(1)


# ----------------------------------------------------------------- Tecplot

def parse_varlocation(header, nvar):
    """Per variable (0-based): True if cell-centred. Unlisted -> NODAL."""
    cell = [False] * nvar
    m = re.search(r"VARLOCATION\s*=\s*\((.*?)\)\s*(,|$)", header, re.I)
    if not m:
        return cell
    for rng, loc in re.findall(r"\[([^\]]*)\]\s*=\s*(\w+)", m.group(1)):
        for part in rng.split(","):
            part = part.strip()
            if "-" in part:
                lo, hi = (int(x) for x in part.split("-"))
            else:
                lo = hi = int(part)
            for v in range(lo, hi + 1):
                if v < 1 or v > nvar:
                    die("VARLOCATION names variable %d of %d" % (v, nvar))
                cell[v - 1] = loc.upper() == "CELLCENTERED"
    return cell


def read_tec(path):
    """A single-zone structured BLOCK file -> dict with the header lines, the
    zone line, the dimensions, the per-variable location and the value
    strings of each variable."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    zones = [n for n, line in enumerate(lines) if line.lstrip().upper().startswith("ZONE")]
    if len(zones) != 1:
        die("%s has %d zones: the splitter takes a single-block case" % (path, len(zones)))
    z = zones[0]
    head = lines[:z]
    vhead = " ".join(head)
    if "VARIABLES" not in vhead.upper():
        die("%s: no VARIABLES line before the zone" % path)
    names = [t.strip('"') for t in re.split(r"[\s,]+", vhead.split("=", 1)[1]) if t.strip('"')]
    zone = lines[z]

    def first_is_number(line):
        tok = line.replace(",", " ").split()
        if not tok:
            return False
        try:
            float(tok[0])
            return True
        except ValueError:
            return False

    if z + 1 < len(lines) and lines[z + 1].strip() and not first_is_number(lines[z + 1]):
        die("%s: the zone header spans several lines (not supported)" % path)
    dims = []
    for key in ("I", "J", "K"):
        m = re.findall(r"\b%s\s*=\s*(\d+)" % key, zone)
        if len(m) != 1:
            die("%s: zone header has %d '%s=' entries: %s" % (path, len(m), key, zone))
        dims.append(int(m[0]))
    I, J, K = dims
    if K < 2 or J < 2 or I < 2:
        die("%s: I=%d J=%d K=%d, every direction needs at least one cell" % (path, I, J, K))
    pack = re.search(r"DATAPACKING\s*=\s*(\w+)", zone, re.I)
    if not pack or pack.group(1).upper() != "BLOCK":
        die("%s: only DATAPACKING=BLOCK is supported" % path)
    cell = parse_varlocation(zone, len(names))

    tokens = []
    for line in lines[z + 1:]:
        tokens.extend(line.replace(",", " ").split())
    n_node, n_cell = I * J * K, (I - 1) * (J - 1) * (K - 1)
    sizes = [n_cell if c else n_node for c in cell]
    if len(tokens) != sum(sizes):
        die("%s: %d values, the header and VARLOCATION need %d" % (path, len(tokens), sum(sizes)))
    for t in tokens:
        try:
            float(t)
        except ValueError:
            die("%s: '%s' is not a number" % (path, t))
    data, off = [], 0
    for n in sizes:
        data.append(tokens[off:off + n])
        off += n
    return dict(head=head, zone=zone, names=names, dims=(I, J, K), cell=cell, data=data)


def zone_header(src, title, I, J, K):
    """The source zone line with only the title and I, J, K replaced."""
    def put_title(m):
        quoted = m.group(2).startswith('"')
        return m.group(1) + ('"%s"' % title if quoted else title)
    out, n = re.subn(r'^(\s*ZONE\s+T\s*=\s*)("[^"]*"|[^,]*)', put_title, src, count=1, flags=re.I)
    if n != 1:
        die("cannot find the zone title (ZONE T = ...) in: %s" % src)
    for key, val in (("I", I), ("J", J), ("K", K)):
        out, n = re.subn(r"\b(%s\s*=\s*)\d+" % key, lambda m, v=val: m.group(1) + str(v), out)
        if n != 1:
            die("cannot set %s in the zone header: %s" % (key, src))
    return out


def split_tec(tec, nbi, nbj, title_fmt):
    """Lines of the split file, and the reassembly check."""
    I, J, K = tec["dims"]
    ni, nj, nk = I - 1, J - 1, K - 1
    si, sj = ni // nbi, nj // nbj
    out = list(tec["head"])
    rebuilt = [[None] * len(v) for v in tec["data"]]
    for bj in range(nbj):
        for bi in range(nbi):
            b = 1 + bi + nbi * bj
            out.append(zone_header(tec["zone"], title_fmt % b, si + 1, sj + 1, K))
            for v, values in enumerate(tec["data"]):
                if tec["cell"][v]:
                    for k in range(nk):
                        for j in range(bj * sj, (bj + 1) * sj):
                            for i in range(bi * si, (bi + 1) * si):
                                g = i + ni * (j + nj * k)
                                out.append(" " + values[g])
                                rebuilt[v][g] = values[g]
                else:
                    for k in range(K):
                        for j in range(bj * sj, (bj + 1) * sj + 1):
                            for i in range(bi * si, (bi + 1) * si + 1):
                                g = i + I * (j + J * k)
                                out.append(" " + values[g])
                                if rebuilt[v][g] is not None and rebuilt[v][g] != values[g]:
                                    die("internal: shared node %d written twice differently" % g)
                                rebuilt[v][g] = values[g]
    # Self-check: the blocks cover the source exactly once (shared nodes aside).
    if rebuilt != tec["data"]:
        die("internal: the blocks do not reassemble into the source field")
    return out


# --------------------------------------------------------------------- bc.txt

def read_bc(path, ni, nj, nk):
    """Records (b, i, j, k, f, type, second-line-or-None), validated."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    recs, n = [], 0
    while n < len(lines):
        if not lines[n].strip():
            if any(line.strip() for line in lines[n:]):
                die("%s: blank line %d inside the records" % (path, n + 1))
            break
        tok = lines[n].replace(",", " ").split()
        if len(tok) != 6:
            die("%s line %d: a record starts with 'b i j k f type': %s" % (path, n + 1, lines[n]))
        b, i, j, k, f, ty = (int(t) for t in tok)
        if ty in CONNECTIVITY:
            die("%s line %d: type %d (%s) carries connectivity the splitter cannot re-map"
                % (path, n + 1, ty, CONNECTIVITY[ty]))
        second = None
        if ty in ONE_VALUE_LINE:
            if n + 1 >= len(lines) or not lines[n + 1].strip():
                die("%s line %d: type %d needs a value line" % (path, n + 1, ty))
            second = lines[n + 1]
            n += 2
        else:
            n += 1
        recs.append((b, i, j, k, f, ty, second))

    # Coverage: every boundary cell of every face exactly once, faces 1..6 in order.
    normal = {1: (0, 1), 2: (0, ni), 3: (1, 1), 4: (1, nj), 5: (2, 1), 6: (2, nk)}
    size = (ni, nj, nk)
    seen, last_f = set(), 0
    for r in recs:
        b, i, j, k, f, ty, _ = r
        if b != 1:
            die("%s: record on block %d, the source must be single-block" % (path, b))
        if f < 1 or f > 6:
            die("%s: face %d" % (path, f))
        if f < last_f:
            die("%s: face %d after face %d, records must be grouped by face in order 1..6" % (path, f, last_f))
        last_f = f
        ijk = (i, j, k)
        ax, val = normal[f]
        if ijk[ax] != val:
            die("%s: record (%d,%d,%d) is not on face %d" % (path, i, j, k, f))
        for a in range(3):
            if not 1 <= ijk[a] <= size[a]:
                die("%s: record (%d,%d,%d) outside the block" % (path, i, j, k))
        if (f, ijk) in seen:
            die("%s: face %d cell (%d,%d,%d) appears twice" % (path, f, i, j, k))
        seen.add((f, ijk))
    expect = 2 * nj * nk + 2 * ni * nk + 2 * ni * nj
    if len(recs) != expect:
        die("%s: %d records, a %dx%dx%d block needs %d" % (path, len(recs), ni, nj, nk, expect))
    return recs


def split_bc(recs, ni, nj, nk, nbi, nbj):
    si, sj = ni // nbi, nj // nbj
    by_face = {f: [r for r in recs if r[4] == f] for f in range(1, 7)}
    out, table = [], {}

    def first(b, i, j, k, f, ty):
        return "%8d%8d%8d%8d%8d%8d" % (b, i, j, k, f, ty)

    def conn(bs, i_s, js, ks, fs):
        return "%8d%8d%8d%8d%8d%8d%8d%8d%8d" % (bs, i_s, js, ks, fs, 1, 0, 0, 1)

    for bj in range(nbj):
        for bi in range(nbi):
            b = 1 + bi + nbi * bj
            i0, j0 = bi * si, bj * sj
            count = 0
            for f in range(1, 7):
                for (_, i, j, k, _, ty, second) in by_face[f]:
                    # tangential range of this block (the normal index is the face's own)
                    if f in (3, 4, 5, 6) and not i0 < i <= i0 + si:
                        continue
                    if f in (1, 2, 5, 6) and not j0 < j <= j0 + sj:
                        continue
                    li = 1 if f == 1 else (si if f == 2 else i - i0)
                    lj = 1 if f == 3 else (sj if f == 4 else j - j0)
                    lk = k
                    outer = ((f == 1 and bi == 0) or (f == 2 and bi == nbi - 1) or
                             (f == 3 and bj == 0) or (f == 4 and bj == nbj - 1) or f in (5, 6))
                    if outer:
                        out.append(first(b, li, lj, lk, f, ty))
                        if second is not None:
                            out.append(second)
                        table[(b, li, lj, lk, f)] = (ty, None)
                    else:
                        if f == 1:
                            tgt = (b - 1, si, lj, lk, 2)
                        elif f == 2:
                            tgt = (b + 1, 1, lj, lk, 1)
                        elif f == 3:
                            tgt = (b - nbi, li, sj, lk, 4)
                        else:
                            tgt = (b + nbi, li, 1, lk, 3)
                        out.append(first(b, li, lj, lk, f, 101))
                        out.append(conn(*tgt))
                        table[(b, li, lj, lk, f)] = (101, tgt)
                    count += 1
            expect = 2 * sj * nk + 2 * si * nk + 2 * si * sj
            if count != expect:
                die("internal: block %d has %d records, %d expected" % (b, count, expect))

    # Self-check: every 101 points at a 101 that points back.
    for key, (ty, tgt) in table.items():
        if ty != 101:
            continue
        back = table.get(tgt)
        if back is None or back[0] != 101 or back[1] != key:
            die("internal: 101 record %s -> %s is not paired" % (key, tgt))
    return out


# ----------------------------------------------------------------- input.ini

def probe_names(sections, path):
    """The probe section names [FUSS-Probes] lists (probe1 = <name>, ...)."""
    names = []
    for name, body in sections:
        if name != "FUSS-Probes":
            continue
        for line in body[1:]:
            s = line.split(";")[0].split("!")[0].split("#")[0]
            m = re.match(r"^\s*probe\d+\s*=\s*(\S+)\s*$", s)
            if m:
                names.append(m.group(1))
    return names


def reindex_probe(body, path, si, sj, nbi):
    """A probe section with its index-position "1 i j k" moved to the split block."""
    out = []
    for line in body:
        m = re.match(r"^(\s*index-position\s*=\s*)([^;!#]*)(.*)$", line)
        if not m:
            out.append(line)
            continue
        tok = m.group(2).split()
        if len(tok) != 4:
            die("%s: index-position needs 'b i j k': %s" % (path, line))
        b, i, j, k = (int(t) for t in tok)
        if b != 1:
            die("%s: index-position on block %d, the source must be single-block" % (path, b))
        bi, bj = (i - 1) // si, (j - 1) // sj
        new = "%d %d %d %d" % (1 + bi + nbi * bj, i - bi * si, j - bj * sj, k)
        out.append(m.group(1) + new + (" " + m.group(3).lstrip() if m.group(3).strip() else ""))
    return out


def fuss_sections(text, path, si=1, sj=1, nbi=1):
    """The [FUSS-*] sections of an INI file, verbatim, as one string each, and the
    section of every probe [FUSS-Probes] names, its index-position re-indexed."""
    sections, cur = [], None
    for line in text.split("\n"):
        m = re.match(r"^\s*\[([^\]]+)\]", line)
        if m:
            cur = [m.group(1).strip(), [line]]
            sections.append(cur)
        elif cur is not None:
            cur[1].append(line)
    probes = probe_names(sections, path)
    kept = []
    for name, body in sections:
        if name in ("FUSS-Multigrid", "HYDRA-Multigrid"):
            die("%s: [%s] is not supported (coarse-level bc files)" % (path, name))
        if name in probes:
            body = reindex_probe(body, path, si, sj, nbi)
        elif not name.startswith("FUSS-"):
            continue
        # trailing blank and comment-only lines separate this section from the next
        while body and (not body[-1].strip() or body[-1].lstrip()[0] in COMMENT_CHARS):
            body.pop()
        for line in body:
            s = line.lstrip()
            if s and s[0] in COMMENT_CHARS and "=" in s:
                die("%s: comment line with '=' in [%s]: %s" % (path, name, line))
        kept.append("\n".join(body))
    if not kept:
        die("%s: no [FUSS-*] section" % path)
    return kept


# ----------------------------------------------------------------------- main

def label(path):
    """A path relative to the repository when it lies inside it."""
    rel = os.path.relpath(os.path.realpath(path), os.path.realpath(ROOT))
    return path if rel.startswith("..") else rel


def write_lines(path, lines):
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def apply_sets(sections, sets):
    """Replace KEY's value in section [SECTION] for every (SECTION, KEY, VALUE)."""
    out = list(sections)
    for sec, key, val in sets:
        hit = [n for n, s in enumerate(out) if re.match(r"^\s*\[%s\]" % re.escape(sec), s)]
        if len(hit) != 1:
            die("--set: %d kept sections [%s]" % (len(hit), sec))
        pat = re.compile(r"^(\s*%s\s*=\s*)([^;!#\n]*?)(\s*([;!#].*)?)$" % re.escape(key), re.M)
        new, n = pat.subn(lambda m: m.group(1) + val + m.group(3), out[hit[0]])
        if n != 1:
            die("--set: [%s] has %d '%s' lines" % (sec, n, key))
        out[hit[0]] = new
    return out


def main(argv):
    sets, rest = [], [argv[0]]
    a = 1
    while a < len(argv):
        if argv[a] == "--set":
            if a + 1 >= len(argv):
                die("--set needs SECTION:KEY=VALUE")
            m = re.match(r"^([^:=]+):([^=]+)=(.+)$", argv[a + 1])
            if not m:
                die("--set %s: expected SECTION:KEY=VALUE" % argv[a + 1])
            sets.append((m.group(1).strip(), m.group(2).strip(), m.group(3).strip()))
            a += 2
        else:
            rest.append(argv[a])
            a += 1
    argv = rest
    if len(argv) != 5:
        sys.stderr.write(__doc__.split("WHY")[0])
        return 2
    src, out = argv[1], argv[4]
    try:
        nbi, nbj = int(argv[2]), int(argv[3])
    except ValueError:
        die("nbi and nbj are integers")
    if nbi < 1 or nbj < 1 or nbi * nbj < 2:
        die("nbi x nbj = %d x %d: at least two blocks" % (nbi, nbj))
    if not os.path.isdir(src):
        die("no case directory %s" % src)
    if os.path.realpath(src) == os.path.realpath(out):
        die("out_dir is the source case")

    with open(os.path.join(src, "input.ini")) as fh:
        env = coupled_env_file(fh.read(), src)
    if env is not None:
        COPY_EXTRA.add(env)
    extra = sorted(f for f in set(os.listdir(os.path.join(src, "INPUT"))) - set(KNOWN_INPUT)
                   if not copied_input(f))
    if extra:
        die("INPUT/ holds files the splitter cannot split or does not know: %s" % " ".join(extra))
    # case.json and expected.json are a generator's description and a verify.py's budgets, not solver inputs
    stray = sorted(f for f in os.listdir(src) if f.endswith((".dat", ".json"))
                   and f not in ("case.json", "expected.json"))
    if stray:
        die("case-root time/probe files are not supported: %s" % " ".join(stray))

    ic = read_tec(os.path.join(src, "INPUT", "ic.tec"))
    if [n.lower() for n in ic["names"][:3]] != ["x", "y", "z"] or any(ic["cell"][:3]) \
            or not all(ic["cell"][3:]):
        die("ic.tec: expected nodal x y z, then cell-centred variables")
    I, J, K = ic["dims"]
    ni, nj, nk = I - 1, J - 1, K - 1
    if ni % nbi or nj % nbj:
        die("%d x %d cells do not divide into %d x %d equal blocks" % (ni, nj, nbi, nbj))

    mesh_path = os.path.join(src, "MESH", "mesh.tec")
    mesh = None
    if os.path.exists(mesh_path):
        mesh = read_tec(mesh_path)
        if mesh["dims"] != ic["dims"]:
            die("MESH/mesh.tec is %s, INPUT/ic.tec %s" % (mesh["dims"], ic["dims"]))
        if [[float(t) for t in v] for v in mesh["data"][:3]] != \
                [[float(t) for t in v] for v in ic["data"][:3]]:
            sys.stderr.write("split_blocks.py: NOTE: MESH/mesh.tec and INPUT/ic.tec hold different "
                             "node coordinates; each is split as it is\n")

    recs = read_bc(os.path.join(src, "INPUT", "bc.txt"), ni, nj, nk)

    with open(os.path.join(src, "input.ini")) as fh:
        sections = fuss_sections(fh.read(), os.path.join(src, "input.ini"), ni // nbi, nj // nbj, nbi)
    with open(os.path.join(src, "FUSS.sh")) as fh:
        sh = fh.read()
    if len(re.findall(r"^MASTERDIR=", sh, re.M)) != 1:
        die("FUSS.sh: expected exactly one MASTERDIR= line")

    # ---- everything read and checked: write
    os.makedirs(os.path.join(out, "INPUT"), exist_ok=True)
    write_lines(os.path.join(out, "INPUT", "ic.tec"), split_tec(ic, nbi, nbj, "B%d-SP"))
    if mesh is not None:
        os.makedirs(os.path.join(out, "MESH"), exist_ok=True)
        write_lines(os.path.join(out, "MESH", "mesh.tec"), split_tec(mesh, nbi, nbj, "Block%d"))
    write_lines(os.path.join(out, "INPUT", "bc.txt"), split_bc(recs, ni, nj, nk, nbi, nbj))
    for name in sorted(os.listdir(os.path.join(src, "INPUT"))):
        if copied_input(name):
            shutil.copyfile(os.path.join(src, "INPUT", name), os.path.join(out, "INPUT", name))

    # MASTERDIR: the repository relative to out_dir (../../.. for test/moving-mesh/<case>)
    master = os.path.relpath(os.path.realpath(ROOT), os.path.realpath(out))
    if not all(p == ".." for p in master.split(os.sep)):
        master = os.path.realpath(ROOT)          # out_dir is not below the repository
        sys.stderr.write("split_blocks.py: NOTE: %s is not below the repository; "
                         "MASTERDIR set to the absolute path %s\n" % (out, master))
    sh = re.sub(r"^MASTERDIR=.*$", "MASTERDIR=" + master, sh, count=1, flags=re.M)
    dst = os.path.join(out, "FUSS.sh")
    with open(dst, "w") as fh:
        fh.write(sh)
    os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    sections = apply_sets(sections, sets)
    header = ("; Written by test/tools/split_blocks.py from %s, cut into %d x %d blocks"
              " (only its FUSS sections are kept)." % (label(src), nbi, nbj))
    for sec, key, val in sets:
        header += "\n; set by --set: [%s] %s %s" % (sec, key, val)
    if "=" in header:
        die("internal: the header comment holds '=' (%s)" % header)
    with open(os.path.join(out, "input.ini"), "w") as fh:
        fh.write(header + "\n\n" + "\n\n".join(sections) + "\n")

    si, sj = ni // nbi, nj // nbj
    nconn = 2 * nk * ((nbi - 1) * nj + (nbj - 1) * ni)
    print("%s -> %s: %d x %d blocks of %d x %d x %d cells, %d boundary records (%d of type 101)"
          % (label(src), label(out), nbi, nbj, si, sj, nk,
             nbi * nbj * (2 * sj * nk + 2 * si * nk + 2 * si * sj), nconn))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
