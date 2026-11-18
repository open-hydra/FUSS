#!/usr/bin/env python3
"""Static check of the MPI tag table (plan 14 S5, D4): every point-to-point call in src/ (MPI_SEND, MPI_RECV,
MPI_ISEND, MPI_IRECV, MPI_SEND_INIT, MPI_RECV_INIT) must carry, as its tag, a named constant TAG_* -- directly, or
through a local variable / dummy argument whose every assignment / actual argument is one; no arithmetic (no rank,
no block number). Every TAG_* constant declared with `parameter` must have one value, the values must be distinct,
in 0..32767 (the least MPI_TAG_UB the standard allows), and below 7001 (hydra-MF's own tags 7001-7008 share
MPI_COMM_WORLD). A tag named TAG_* that is not a declared integer parameter, a parameter whose value is an expression
(TAG_X = 7000 + 5), and any point-to-point call this check does not resolve (SENDRECV, SSEND, BSEND, RSEND, ISSEND,
PROBE, IPROBE, MPROBE, MRECV) are violations too.
Usage: python3 -I tag_table_check.py <repo> <tree-ish>     (git mode)
       python3 -I tag_table_check.py --dir <src directory>  (files on disk, e.g. the stage worktree's src/)
Prints the table and every call with its resolved constant; exit 1 on any violation."""
import re, sys, subprocess, os
P2P = re.compile(r'\bcall\s+MPI_(I?SEND|I?RECV|SEND_INIT|RECV_INIT)\s*\(', re.I)
OTHER = re.compile(r'\bcall\s+MPI_(SENDRECV\w*|I?[BSR]SEND\w*|ISSEND|I?M?PROBE|I?MRECV)\s*\(', re.I)   # not resolved here: any is a violation

def load():
    out = {}
    if sys.argv[1] == '--dir':
        for root, _, names in os.walk(sys.argv[2]):
            for n in names:
                if n.endswith('.f90'):
                    p = os.path.join(root, n); out[os.path.relpath(p, sys.argv[2])] = open(p, errors='replace').read()
    else:
        repo, pin = sys.argv[1], sys.argv[2]
        git = lambda *a: subprocess.run(['git', '-C', repo] + list(a), stdout=subprocess.PIPE,
                                        universal_newlines=True, encoding='utf-8', check=True).stdout
        for f in git('ls-tree', '-r', '--name-only', pin, 'src').split('\n'):
            if f.endswith('.f90'):
                out[f] = git('show', f'{pin}:{f}')
    return out

def logical_lines(txt):
    """(first physical line number, joined statement) with comments stripped and & continuations joined."""
    res, buf, start = [], '', None
    for n, raw in enumerate(txt.split('\n'), 1):
        line = re.sub(r'!.*$', '', raw) if not re.match(r'\s*#', raw) else ''
        if not line.strip() and not buf: continue
        if start is None: start = n
        s = line.rstrip()
        cont = s.endswith('&')
        s = s[:-1] if cont else s
        s = re.sub(r'^\s*&', '', s) if buf else s
        buf += ' ' + s.strip()
        if not cont:
            for part in buf.split(';'):          # statements joined by ';' (e.g. case (X); h%tag = TAG_X)
                if part.strip(): res.append((start, part.strip()))
            buf, start = '', None
    return res

def split_args(s):
    depth, cur, args = 0, '', []
    for ch in s:
        if ch == '(' : depth += 1
        if ch == ')' :
            if depth == 0: break
            depth -= 1
        if ch == ',' and depth == 0: args.append(cur.strip()); cur = ''
        else: cur += ch
    args.append(cur.strip()); return args

files = load()
consts, errors = {}, []
units = []   # (file, unit name, dummy list, [(line, stmt)])
for f, txt in files.items():
    cur = None
    for n, st in logical_lines(txt):
        if re.search(r'\bparameter\b', st, re.I):
            for name, val in re.findall(r'\b(TAG_[A-Z0-9_]+)\s*=\s*(-?\d+)\s*(?=,|$)', st):   # a literal only; widen deliberately if upstream writes 5_I4
                consts.setdefault(name, set()).add((int(val), f, n))
        m = re.match(r'(?:(?:recursive|pure|elemental|module)\s+)*(subroutine|function)\s+(\w+)\s*(\(([^)]*)\))?', st, re.I)
        if m:
            cur = [f, m.group(2).lower(), [a.strip().lower() for a in (m.group(4) or '').split(',') if a.strip()], []]
            units.append(cur); continue
        if cur is not None: cur[3].append((n, st))

values = {}
for name, s in sorted(consts.items()):
    vals = {v for v, _, _ in s}
    if len(vals) > 1: errors.append(f'{name} declared with several values: {sorted(s)}')
    values[name] = min(vals)
    if not (0 <= values[name] <= 32767): errors.append(f'{name} = {values[name]} outside 0..32767')
    if values[name] >= 7001: errors.append(f'{name} = {values[name]} not below 7001 (hydra-MF tags)')
seen = {}
for name, v in values.items():
    if v in seen: errors.append(f'{name} and {seen[v]} share the value {v}')
    seen[v] = name

def resolve(expr, unit, depth=0):
    """set of TAG_* names expr can take, or None if some path is not a constant."""
    e = expr.strip()
    if re.fullmatch(r'TAG_[A-Z0-9_]+', e): return {e} if e in values else None   # a declared integer parameter only
    if not re.fullmatch(r'[A-Za-z_][\w%]*', e) or depth > 4: return None
    var = e.lower(); got = set(); found = False
    for n, st in unit[3]:
        m = re.match(rf'(?:if\s*\(.*\)\s*)?{re.escape(var)}\s*=\s*(.+)$', st, re.I)
        if m and not re.match(r'\s*=', m.group(1)):
            found = True
            r = resolve(m.group(1), unit, depth + 1)
            if r is None: return None
            got |= r
    if var in unit[2]:   # dummy argument: every call site's actual argument at that position
        pos = unit[2].index(var); found = True
        for u in units:
            for n, st in u[3]:
                m = re.search(rf'\bcall\s+{unit[1]}\s*\((.*)\)', st, re.I)
                if not m: continue
                a = split_args(m.group(1) + ')')
                if pos >= len(a): return None
                r = resolve(a[pos], u, depth + 1)
                if r is None: return None
                got |= r
    return got if found and got else None

calls = 0
for u in units:
    for n, st in u[3]:
        if OTHER.search(st): errors.append(f'{u[0]}:{n} [{u[1]}] point-to-point call this check does not resolve: {st[:60]}')
        m = P2P.search(st)
        if not m: continue
        calls += 1
        args = split_args(st[m.end():])
        tag = args[4] if len(args) > 4 else '?'
        r = resolve(tag, u)
        mark = 'ok ' if r else 'BAD'
        if not r: errors.append(f'{u[0]}:{n} [{u[1]}] tag "{tag}" is not a named TAG_* constant')
        print(f'{mark} {u[0]}:{n} [{u[1]}] MPI_{m.group(1).upper()} tag {tag} -> {sorted(r) if r else "-"}')
if calls == 0: errors.append('no point-to-point MPI call found: wrong directory or tree-ish')
print(f'\n{calls} point-to-point calls; tag constants:')
for name in sorted(values, key=values.get):
    print(f'  {values[name]:6d}  {name}  ({", ".join(sorted({f for _, f, _ in consts[name]}))})')
if errors:
    print('\nVIOLATIONS:'); [print('  ' + e) for e in errors]; sys.exit(1)
print('\nPASS: every tag a named constant, values distinct, all in 0..7000')
