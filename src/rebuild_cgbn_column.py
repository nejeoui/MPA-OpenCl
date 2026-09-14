#!/usr/bin/env python3
"""Rebuild a report's CGBN column, and its head-to-head section, from the CSV.

Two things get rebuilt. The CGBN column comes from a cgbn_results TSV, for a
report written against the wrong file. The whole of section 5 is then recomputed
from the CSV's per-row ops/s, because GPU_Host wrote that section without
normalising the GPU time to the CPU item count - so every figure there was out
by the ratio of the two counts whenever they differed. Section 4 and the CSV
were always correct and are left alone.
"""
import csv, re, sys, os

def rate(ops):
    if ops is None or ops <= 0: return "n/a"
    if ops >= 1e9: return f"{ops/1e9:.2f} G"
    if ops >= 1e6: return f"{ops/1e6:.2f} M"
    if ops >= 1e3: return f"{ops/1e3:.2f} k"
    return f"{ops:.1f}"

def load_tsv(path):
    dev, d = None, {}
    for line in open(path):
        if line.startswith('#'):
            m = re.search(r'device (.*?), items=(\d+)', line)
            if m: dev = m.group(1)
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 4:
            d[(f[0], f[1])] = (int(f[2]), float(f[3]))
    return dev, d

def main(report_base, tsv_path):
    dev, cg = load_tsv(tsv_path)
    csv_path, md_path = report_base + '_Report.csv', report_base + '_Report.md'

    rows = list(csv.reader(open(csv_path)))
    header = rows[0]
    devices = {r[1] for r in rows[1:] if r and r[0].startswith('opencl')}
    if dev not in devices:
        sys.exit(f"refusing: TSV was produced on '{dev}', report devices are {devices}")

    # best MPA ops/s per (modulus, operation), for the ratio column
    best = {}
    for r in rows[1:]:
        if r[0] == 'opencl-kernel' and not r[1].startswith('cpu') and float(r[8]) > 0:
            k = (r[4], r[6]); v = float(r[9])
            if v > best.get(k, 0): best[k] = v

    n = 0
    for r in rows[1:]:
        if r[0] == 'library' and r[3] == 'cgbn':
            k = (r[4], r[6])
            if k in cg:
                items, secs = cg[k]
                r[7], r[8], r[9] = str(items), f"{secs:.9f}", f"{items/secs:.3f}"
                n += 1
    with open(csv_path, 'w', newline='') as fh:
        csv.writer(fh).writerows(rows)
    print(f"{csv_path}: {n} cgbn rows rewritten")

    md = open(md_path).read()
    md = re.sub(r'\| CGBN \| .*? \|',
                f'| CGBN | {len(cg)} rows from `{os.path.basename(tsv_path)}` |', md, count=1)

    # Section 5 is rebuilt entirely from the CSV: per (modulus, operation) the
    # best GPU and best CPU-OpenCL variant, the three CPU libraries, CGBN, and
    # the four ratios - all from ops/s, which carries each row's own item count.
    gbest, cbest, lib = {}, {}, {}
    for r in rows[1:]:
        if not r: continue
        k = (r[4], r[6])
        if r[0] == 'opencl-kernel' and float(r[8]) > 0:
            tgt = cbest if r[1].startswith('cpu') else gbest
            v = float(r[9])
            if v > tgt.get(k, (0, ''))[0]: tgt[k] = (v, r[3])
        elif r[0] == 'library' and float(r[8]) > 0:
            lib[(r[3],) + k] = float(r[9])

    out, modulus, section = [], None, None
    for line in md.split('\n'):
        h = re.match(r'#{3,4} (.+?) \(\d+-bit\)', line)
        if h: modulus = h.group(1)
        if line.startswith('## 4.'): section = 4
        elif line.startswith('## 5.'): section = 5
        elif line.startswith('## 6.'): section = 6
        if line.startswith('| ') and modulus and section in (4, 5) and '---' not in line:
            cells = [c.strip() for c in line.strip().strip('|').split('|')]
            op = cells[0]
            if op != 'Operation' and section == 4 and (modulus, op) in cg:
                items, secs = cg[(modulus, op)]
                cells[-1] = rate(items / secs)
                line = '| ' + ' | '.join(cells) + ' |'
            elif op != 'Operation' and section == 5 and len(cells) == 13:
                k = (modulus, op)
                g, gv = gbest.get(k, (0, 'none'))
                c, cv = cbest.get(k, (0, 'none'))
                cops = 0
                if k in cg:
                    items, secs = cg[k]; cops = items / secs
                rat = lambda a, b: f"{a/b:.2f}x" if a > 0 and b > 0 else "n/a"
                cells = [op, gv if g else 'none', rate(g) if g else 'n/a',
                         cv if c else 'none', rate(c) if c else 'n/a',
                         rate(lib.get(('gmp-1t',) + k, 0)),
                         rate(lib.get(('gmp-nt',) + k, 0)),
                         rate(lib.get(('openssl-nt',) + k, 0)),
                         rate(cops),
                         rat(g, c), rat(g, lib.get(('gmp-nt',) + k, 0)),
                         rat(g, lib.get(('openssl-nt',) + k, 0)), rat(g, cops)]
                line = '| ' + ' | '.join(cells) + ' |'
        out.append(line)
    md = '\n'.join(out)

    banner = (f"> **CGBN column rebuilt offline.** The run that produced this report loaded the\n"
              f"> wrong `cgbn_results` file. Its CGBN column, and the `GPU vs CGBN` ratios, were\n"
              f"> regenerated from `{os.path.basename(tsv_path)}` (device `{dev}`, {len(cg)} rows)\n"
              f"> by `rebuild_cgbn_column.py`, which also recomputed section 5 from the CSV:\n"
              f"> as written, that section did not normalise the GPU time to the CPU item\n"
              f"> count, so its figures were out by the ratio of the two. Section 4, and every\n"
              f"> measured value in the CSV, are untouched.\n\n")
    md = re.sub(r'(^# .*?\n\n)', r'\1' + banner, md, count=1, flags=re.S)
    open(md_path, 'w').write(md)
    print(f"{md_path}: CGBN cells and ratios rebuilt, provenance banner added")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("usage: rebuild_cgbn_column.py <Report base name> <cgbn_results TSV>")
    main(sys.argv[1], sys.argv[2])
