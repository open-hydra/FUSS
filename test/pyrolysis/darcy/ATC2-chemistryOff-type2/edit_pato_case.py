#!/usr/bin/env python3
"""edit_pato_case.py <case-dir> <tacot-noeps-dir>

Plan 06 section 4.5 edits (1)-(3) on a copy of PATO's
AblationTestCase_2.x_chemistryOff.  Reproducible and idempotent;
Python 3.6, standard library only.

(1) constant/porousMat/BoundaryConditions: every data row keeps its
    time token verbatim and gets p_total_w, rhoUeCH, h_r from the row
    at t = 0.1 (tokens taken literally from that row), then
    hconv = 300, Tedge = 2800, chemistryOn = 0.  Tokens are replaced
    in place, so the separators of each row are preserved.
(2) constant/porousMat/porousMatProperties: MaterialPropertiesDirectory
    -> absolute path of <tacot-noeps-dir>.  GasPropertiesFile untouched.
(3) system/controlDict: endTime -> 60.
(4) system/porousMat/plotDict: `p` added after `Ta` to the probed
    `fields` (plan 07 section 4.4: output/porousMat/scalar/p_plot then
    holds the gas pressure at the same seven depths; Ta_plot is
    unchanged byte for byte).
movingMesh is left as shipped.
"""
import os
import re
import sys

TOKEN = re.compile(r'\S+')


def strip_comments(lines):
    """Mirror PATO's IOFunctions.H removeComments: '//' to end of line,
    '/* ... */' possibly spanning lines.  Returns the code part of each line."""
    out = []
    inside = False
    for line in lines:
        line = line.split('//', 1)[0]
        s = line.find('/*')
        e = line.find('*/')
        if s >= 0:
            inside = True
        if inside:
            if s >= 0:
                if e >= 0 and e > s:
                    line = line[:s] + line[e + 2:]
                    inside = False
                else:
                    line = line[:s]
            else:
                if e >= 0:
                    line = line[e + 2:]
                    inside = False
                else:
                    line = ''
        out.append(line)
    return out


def edit_boundary_conditions(path):
    with open(path, 'r', newline='') as f:
        text = f.read()
    lines = text.split('\n')
    code = strip_comments(lines)
    data_idx = []
    for i, c in enumerate(code):
        tok = c.split()
        if not tok:
            continue
        try:
            float(tok[0])
        except ValueError:
            continue
        data_idx.append(i)
    if not data_idx:
        raise SystemExit('no data rows in ' + path)
    # the reference row: t = 0.1
    ref = None
    for i in data_idx:
        tok = code[i].split()
        if abs(float(tok[0]) - 0.1) < 1e-12:
            ref = tok
            break
    if ref is None:
        raise SystemExit('no row at t = 0.1 in ' + path)
    if len(ref) != 7:
        raise SystemExit('expected 7 columns, found %d in %r' % (len(ref), ref))
    p_tok, rhoUeCH_tok, h_r_tok = ref[1], ref[2], ref[3]
    new_tail = [p_tok, rhoUeCH_tok, h_r_tok, '300', '2800', '0']
    for i in data_idx:
        line = lines[i]
        if code[i] != line:
            raise SystemExit('data row with a comment, not handled: %r' % line)
        spans = [m.span() for m in TOKEN.finditer(line)]
        if len(spans) != 7:
            raise SystemExit('expected 7 tokens, found %d in %r' % (len(spans), line))
        # replace tokens 1..6 in place, from the right so spans stay valid
        for k in range(6, 0, -1):
            a, b = spans[k]
            line = line[:a] + new_tail[k - 1] + line[b:]
        lines[i] = line
    new = '\n'.join(lines)
    with open(path, 'w', newline='') as f:
        f.write(new)
    return new != text


def edit_material_dir(path, tacot_dir):
    with open(path, 'r', newline='') as f:
        text = f.read()
    pat = re.compile(r'^(\s*MaterialPropertiesDirectory\s+)"[^"]*"(\s*;)', re.M)
    if len(pat.findall(text)) != 1:
        raise SystemExit('MaterialPropertiesDirectory not found exactly once in ' + path)
    new = pat.sub(lambda m: m.group(1) + '"' + tacot_dir + '"' + m.group(2), text)
    with open(path, 'w', newline='') as f:
        f.write(new)
    return new != text


def edit_end_time(path, end_time='60'):
    with open(path, 'r', newline='') as f:
        text = f.read()
    pat = re.compile(r'^(endTime\s+)\S+(\s*;)', re.M)
    if len(pat.findall(text)) != 1:
        raise SystemExit('endTime not found exactly once in ' + path)
    new = pat.sub(lambda m: m.group(1) + end_time + m.group(2), text)
    with open(path, 'w', newline='') as f:
        f.write(new)
    return new != text


def edit_plot_fields(path):
    with open(path, 'r', newline='') as f:
        text = f.read()
    pat = re.compile(r'^(fields\s*\(\s*)([^)]*)(\))', re.M)
    m = pat.search(text)
    if m is None or len(pat.findall(text)) != 1:
        raise SystemExit('a single `fields ( ... )` list not found in ' + path)
    names = m.group(2).split()
    if 'Ta' not in names:
        raise SystemExit('`Ta` is not a probed field in ' + path)
    if 'p' in names:
        return False
    body = m.group(2)
    i = body.index('Ta') + len('Ta')
    eol = body.find('\n', i)
    ta = m.start(2) + body.index('Ta')                     # the indent of the `Ta` line, from the file
    indent = re.match(r'[ \t]*', text[text.rfind('\n', 0, ta) + 1:]).group(0)
    body = body[:eol + 1] + indent + 'p\n' + body[eol + 1:] if eol >= 0 else body[:i] + ' p' + body[i:]
    new = text[:m.start(2)] + body + text[m.end(2):]
    with open(path, 'w', newline='') as f:
        f.write(new)
    return True


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    case = os.path.abspath(sys.argv[1])
    tacot = os.path.abspath(sys.argv[2])
    for name in ('virgin', 'char', 'constantProperties', 'gasProperties'):
        if not os.path.isfile(os.path.join(tacot, name)):
            raise SystemExit('missing %s in %s' % (name, tacot))
    bc = os.path.join(case, 'constant', 'porousMat', 'BoundaryConditions')
    pm = os.path.join(case, 'constant', 'porousMat', 'porousMatProperties')
    cd = os.path.join(case, 'system', 'controlDict')
    pd = os.path.join(case, 'system', 'porousMat', 'plotDict')
    print('(1) BoundaryConditions   %s' % ('edited' if edit_boundary_conditions(bc) else 'already as required'))
    print('(2) porousMatProperties  %s' % ('edited' if edit_material_dir(pm, tacot) else 'already as required'))
    print('(3) controlDict          %s' % ('edited' if edit_end_time(cd) else 'already as required'))
    print('(4) plotDict             %s' % ('edited' if edit_plot_fields(pd) else 'already as required'))


if __name__ == '__main__':
    main()
