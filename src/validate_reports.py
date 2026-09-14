import csv, glob, os, re, sys, statistics as st

V = ['w8','w16','w32','w32-opt','w32-o64','w32-il','w32-il64']
FULL = {v: (35 if v == 'w32' else 75) for v in V}

def analyse(csv_path):
    md = csv_path[:-4] + '.md'
    txt = open(md, errors='replace').read() if os.path.exists(md) else ''
    def g(pat, d=None):
        m = re.search(pat, txt)
        return m.group(1).strip() if m else d
    r = {'name': os.path.basename(csv_path)[:-11], 'notes': []}
    r['cu']      = g(r'\| Compute units \| (\d+) \|', '?')
    r['runtime'] = g(r'\| OpenCL version \| (.*?) \|', '?')
    r['threads'] = int(g(r'\| OpenMP threads used \| (\d+) \|', '0') or 0)
    r['cpufix']  = 'Cost weighting drives the wide cells' in txt

    cov, mm, inf, cells = {}, 0, 0, 0
    K, C, L = {}, {}, {}
    cgbn_items, sweep_items = set(), set()
    for row in csv.DictReader(open(csv_path)):
        kind, dev = row['kind'], row['device']
        secs = float(row['seconds'])
        key = (row['modulus'], row['operation'])
        if kind == 'opencl-kernel' and not dev.startswith('cpu'):
            cells += 1; mm += int(row['mismatches'])
            if secs <= 0: inf += 1
            else:
                cov[row['kernel']] = cov.get(row['kernel'], 0) + 1
                K.setdefault(key, {})[row['kernel']] = float(row['ops_per_sec'])
                sweep_items.add(row['items'])
        elif kind == 'library':
            if row['kernel'] == 'cgbn':
                C[key] = row['seconds']; cgbn_items.add(row['items'])
            elif secs > 0:
                L[(row['kernel'],) + key] = float(row['ops_per_sec'])
    r['cells'], r['mismatches'], r['inf'] = cells, mm, inf
    r['variants'] = sum(1 for v in V if cov.get(v, 0) == FULL[v])
    r['K'], r['C'] = K, C

    # --- MPA kernel data ---------------------------------------------------
    r['mpa_ok'] = mm == 0 and cells > 0
    if mm:  r['notes'].append(f'{mm} mismatched items')
    if inf: r['notes'].append(f'{inf} zero-second rows (filter seconds>0)')
    if len(sweep_items) and max(int(i) for i in sweep_items) < 20000:
        r['notes'].append('all cells below 20000 items: launch-overhead bound')
        r['mpa_ok'] = False

    # --- CPU library baselines --------------------------------------------
    bad = [f'{m}/{o}' for m in ('p1024','p2048')
           for o in ('REDUCE','DIVIDE','ISQRT','MODMUL','MODEXP')
           if ('gmp-1t',m,o) in L and ('gmp-nt',m,o) in L
           and L[('gmp-nt',m,o)] < L[('gmp-1t',m,o)]]
    r['cpu_ok'] = r['cpufix'] and not bad
    if not r['cpufix']: r['notes'].append('CPU baselines predate the 2026-09-12 timing fix')
    elif bad:           r['notes'].append(f'{len(bad)} CPU cells where nT < 1T ({r["threads"]} threads)')

    # --- CGBN column -------------------------------------------------------
    r['cgbn_ok'] = bool(C)
    if not C:
        r['cgbn_ok'] = False; r['notes'].append('no CGBN column')
    else:
        if cgbn_items != sweep_items:
            r['cgbn_ok'] = False
            r['notes'].append(f'CGBN at {sorted(cgbn_items)} items vs sweep at {sorted(sweep_items)}')
    return r

reports = [analyse(f) for f in sorted(glob.glob('*_Report.csv'))]

# contamination: byte-identical CGBN columns between two reports
for i, a in enumerate(reports):
    for b in reports[i+1:]:
        common = set(a['C']) & set(b['C'])
        if len(common) >= 30 and all(a['C'][k] == b['C'][k] for k in common):
            for x, y in ((a, b), (b, a)):
                x['cgbn_ok'] = False
                x['notes'].append(f'CGBN column identical to {y["name"]} (contaminated)')

print(f"{'report':<32}{'vars':>5}{'MPA':>5}{'CPU':>5}{'CGBN':>6}  notes")
for r in sorted(reports, key=lambda r: (-r['variants'], r['name'])):
    y = lambda b: ' ok ' if b else 'BAD '
    print(f"{r['name']:<32}{r['variants']:>4}/7{y(r['mpa_ok']):>5}{y(r['cpu_ok']):>5}{y(r['cgbn_ok']):>6}  "
          + '; '.join(r['notes']))

clean = [r for r in reports if r['mpa_ok'] and r['cpu_ok'] and r['cgbn_ok'] and r['variants'] == 7]
print(f"\nfully usable (7/7 variants, MPA + CPU + CGBN all valid): "
      + (', '.join(r['name'] for r in clean) or 'none'))
print(f"MPA data usable: {sum(1 for r in reports if r['mpa_ok'])}/{len(reports)} reports")
