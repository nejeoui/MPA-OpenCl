#!/usr/bin/env python3
"""Fold cgbn_results.tsv into an already-generated GPU_Host report.

    ./cgbn_merge.py NVIDIA_GeForce_RTX_4070_Ti_Report.md [cgbn_results.tsv]

Appends a "CGBN comparison" section to the .md and cgbn rows to the .csv,
both idempotently. Use it when CGBN was measured after the sweep and you do
not want to pay for the whole sweep again.
"""
import csv, os, sys

MARK = "## CGBN comparison"
OVERHEAD_S = 20e-6          # below this a CUDA kernel time is mostly launch cost
STARVED_N  = 4096           # below this an item count cannot fill a modern GPU


def rate(ops):
    if ops is None or ops <= 0:
        return "n/a"
    for div, suf in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if ops >= div:
            return f"{ops/div:.2f} {suf}"
    return f"{ops:.1f}"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    md = sys.argv[1]
    tsv = sys.argv[2] if len(sys.argv) > 2 else "cgbn_results.tsv"
    csvp = md[:-3] + ".csv" if md.endswith(".md") else md + ".csv"
    for p in (md, tsv, csvp):
        if not os.path.exists(p):
            sys.exit(f"missing: {p}")

    cgbn, order = {}, []
    with open(tsv) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) != 4:
                continue
            mod, op, items, secs = parts[0], parts[1], int(parts[2]), float(parts[3])
            cgbn[(mod, op)] = (items, secs)
    if not cgbn:
        sys.exit(f"{tsv} holds no data rows")

    best, mods, ops, bits = {}, [], [], {}
    with open(csvp) as f:
        for r in csv.DictReader(f):
            if r["kind"] != "opencl-kernel" or r["device_type"] != "GPU":
                continue
            if int(r["mismatches"]):
                continue
            k = (r["modulus"], r["operation"])
            v = float(r["ops_per_sec"])
            if k not in best or v > best[k][0]:
                best[k] = (v, r["kernel"], int(r["items"]))
            if r["modulus"] not in mods:
                mods.append(r["modulus"]); bits[r["modulus"]] = r["bits"]
            if r["operation"] not in ops:
                ops.append(r["operation"])

    out = [MARK, "",
           f"Merged from `{tsv}` after the sweep, by `cgbn_merge.py`.",
           "MPA columns are the fastest correct kernel on the GPU, kernel time only.",
           "", ]
    flagged = set()
    for m in mods:
        rows = []
        for o in ops:
            if (m, o) not in cgbn:
                continue
            citems, csecs = cgbn[(m, o)]
            cops = citems / csecs if csecs > 0 else 0
            b = best.get((m, o))
            cmark = mmark = ""
            if csecs < OVERHEAD_S:
                cmark = " †"; flagged.add("cgbn")
            if b and b[2] < STARVED_N:
                mmark = " ‡"; flagged.add("mpa")
            ratio = f"{cops/b[0]:.2f}x" if b and cops > 0 else "n/a"
            trust = "" if not (cmark or mmark) else " ⚠"
            rows.append(f"| {o} | {b[1] if b else 'n/a'} | "
                        f"{b[2] if b else 0} | "
                        f"{rate(b[0]) if b else 'n/a'}{mmark} | "
                        f"{rate(cops)}{cmark} | {ratio}{trust} |")
        if not rows:
            continue
        out += [f"### {m} ({bits[m]}-bit)", "",
                "| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |",
                "|---|---|---|---|---|---|"] + rows + [""]

    if flagged:
        out += ["**Only rows without a ⚠ are a fair comparison.**", ""]
    if "cgbn" in flagged:
        out += [f"† CGBN kernel time below {OVERHEAD_S*1e6:.0f} us: that is the CUDA launch "
                "floor, not arithmetic. CGBN's real throughput here is higher than shown, "
                "so these rows understate CGBN.", ""]
    if "mpa" in flagged:
        out += [f"‡ MPA ran fewer than {STARVED_N} items. One work-item per thread means a "
                "modern GPU sits near-idle, so these rows understate MPA, and increasingly "
                "so as the modulus grows. Re-run with `--min-items` to fix.", ""]

    body = open(md).read()
    if MARK in body:
        body = body[:body.index(MARK)].rstrip() + "\n\n"
    else:
        body = body.rstrip() + "\n\n"
    open(md, "w").write(body + "\n".join(out).rstrip() + "\n")

    # Drop any CGBN rows from a previous merge before writing the current ones,
    # otherwise a re-measured run leaves stale numbers behind in the csv.
    keep = [l for l in open(csvp).read().splitlines() if ",cgbn," not in l]
    with open(csvp, "w") as f:
        f.write("\n".join(keep) + "\n")
        for (m, o), (items, secs) in cgbn.items():
            f.write(f"library,CGBN,GPU,cgbn,{m},{bits.get(m,'')},{o},"
                    f"{items},{secs:.9f},{items/secs:.3f},0\n")

    print(f"{md}: {sum(1 for k in cgbn if k[0] in mods)} CGBN rows merged")


if __name__ == "__main__":
    main()
