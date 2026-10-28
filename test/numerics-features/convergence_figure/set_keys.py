"""set_keys.py <input.ini> 'Section:key=value;Section:key=value;Section:key='

Sets keys of a FUSS input.ini in place (FiNeR ini syntax). An existing key is replaced, a
missing one is added right after its section header, a missing section is appended at the end.
An empty value removes the key. Comment and blank lines are kept; a line starting with ';' or
'[' ends a section.
"""
import sys


def main(path, spec):
    lines = open(path).read().split('\n')
    for item in filter(None, (s.strip() for s in spec.split(';'))):
        sec, kv = item.split(':', 1)
        key, val = (s.strip() for s in kv.split('=', 1))
        hdr = f'[{sec}]'
        stripped = [l.strip() for l in lines]
        if hdr not in stripped:
            if val:
                lines += ['', hdr, f'{key} = {val}']
            continue
        h = stripped.index(hdr)
        end = next((i for i in range(h + 1, len(lines)) if stripped[i].startswith(('[', ';'))), len(lines))
        for i in range(h + 1, end):
            if lines[i].split('=')[0].strip() == key:
                if val:
                    lines[i] = f'{key} = {val}'
                else:
                    del lines[i]
                break
        else:
            if val:
                lines.insert(h + 1, f'{key} = {val}')
    if lines[-1].strip():
        lines.append('')
    open(path, 'w').write('\n'.join(lines))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
