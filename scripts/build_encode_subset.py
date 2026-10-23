"""Generate the bundled ENCODE H3K4me3 demo signal bigWig.

Takes the downloaded ENCODE file ``ENCFF893ZGR`` (H3K4me3 fold change over
control, mouse mm10, experiment ENCSR000CDQ) and produces a compact
single-window subset for distribution inside ``pychipseeker``:

1. verifies the source file's MD5 against the ENCODE-published value;
2. picks the densest chromosome (highest mean signal among the largest
   chromosomes) and the best ~20 Mb window within it;
3. re-aggregates the signal into ``--bin-bp`` mean bins and drops bins
   below ``--min-value`` (background) so the bundled file stays tiny
   (gaps encode "no data", which is faithful for fold-change tracks);
4. writes the bigWig with ``build_bigwig_py`` (pure Python; the MinGW
   pyBigWig/libBigWig writer crashes on Windows) plus a provenance JSON.

Usage::

    python build_encode_subset.py --source ENCFF893ZGR.bigWig \
        --out-dir ../src/pychipseeker/data
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import urllib.request
from urllib.parse import urlencode

import pyBigWig

SOURCE_MD5 = "f07c99ceead499e3e9d5c4befe92d6e1"
SOURCE_ACCESSION = "ENCFF893ZGR"
SOURCE_SIZE_BYTES = 41_747_495
EXPERIMENT = "ENCSR000CDQ"
ASSEMBLY = "mm10"


def _md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _experiment_metadata():
    """Best-effort fetch of ENCODE experiment metadata (offline-safe)."""
    base = "https://www.encodeproject.org"
    try:
        url = base + f"/experiments/{EXPERIMENT}/?{urlencode({'format': 'json'})}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            exp = json.load(resp)
        return {
            "description": exp.get("description"),
            "biosample": exp.get("biosample_summary"),
            "target": (exp.get("target") or {}).get("title")
            if isinstance(exp.get("target"), dict)
            else exp.get("target"),
        }
    except Exception:  # noqa: BLE001 - metadata is nice-to-have only
        return {
            "description": "H3K4me3 ChIP-seq on 8-week mouse placenta "
            "(fold change over control)",
            "biosample": "Mus musculus strain B6NCrl placenta tissue "
            "female adult (2 months)",
            "target": "H3K4me3 (Mus musculus)",
        }


def _choose_chromosome(bw, chroms):
    """Return the chromosome with the highest mean signal (coarse scan)."""
    candidates = sorted(chroms, key=chroms.get, reverse=True)[:8]
    scores = {}
    for chrom in candidates:
        stats = bw.stats(chrom, 0, chroms[chrom], nBins=200, type="mean")
        valid = [v for v in stats if v is not None]
        scores[chrom] = sum(valid) / len(valid) if valid else 0.0
        print(f"  {chrom}: mean={scores[chrom]:.3f}", flush=True)
    chrom = max(scores, key=scores.get)
    return chrom, scores


def _choose_window(bw, chrom, chrom_size, window):
    """Center a `window` span on the highest local mean-signal region."""
    n_bins = 2000
    step = chrom_size / n_bins
    values = bw.stats(chrom, 0, chrom_size, nBins=n_bins, type="mean")
    half = int(window / step / 2)
    best_start, best_score = 0, -1.0
    for i in range(n_bins):
        lo, hi = max(0, i - half), min(n_bins, i + half)
        segment = [v for v in values[lo:hi] if v is not None]
        if len(segment) < (hi - lo) * 0.6:
            continue
        score = sum(segment) / len(segment) * len(segment) / (hi - lo)
        if score > best_score:
            best_score, best_start = score, int(lo * step)
    return best_start, best_score


def _extract_entries(bw, chrom, start, end, bin_bp, min_value):
    """Yield (chrom, bin_start, bin_end, value) for signal-carrying bins."""
    entries = []
    chunk = 256 << 10  # 256 kb per stats() call
    for lo in range(start, end, chunk):
        hi = min(lo + chunk, end)
        n_bins = (hi - lo) // bin_bp
        if n_bins <= 0:
            continue
        values = bw.stats(chrom, lo, hi, nBins=n_bins, type="mean")
        for i, value in enumerate(values):
            if value is None or value < min_value:
                continue
            bin_start = lo + i * bin_bp
            entries.append((chrom, bin_start, bin_start + bin_bp, round(value, 4)))
    return entries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        default=os.path.join(os.environ.get("TEMP", "."), "ENCFF893ZGR.bigWig"),
        help="path to the downloaded ENCODE bigWig",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="directory for the bundled .bw and .json (pychipseeker/data)",
    )
    parser.add_argument("--window-mb", type=int, default=20, help="window span in Mb")
    parser.add_argument("--bin-bp", type=int, default=25, help="bin size in bp")
    parser.add_argument(
        "--min-value",
        type=float,
        default=0.05,
        help="drop bins whose mean fold change is below this (default 0.05)",
    )
    parser.add_argument(
        "--compress-level", type=int, default=9, help="zlib level (default 9)"
    )
    args = parser.parse_args(argv)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    import sys

    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    from build_bigwig_py import write_bigwig

    print(f"source: {args.source}")
    digest = _md5(args.source)
    if digest != SOURCE_MD5:
        raise SystemExit(
            f"source md5 mismatch: got {digest}, expected {SOURCE_MD5}"
        )
    print("md5 verified:", digest)

    bw = pyBigWig.open(args.source)
    chroms = dict(bw.chroms())
    print(f"chromosomes: {len(chroms)}")
    chrom, _scores = _choose_chromosome(bw, chroms)
    chrom_size = chroms[chrom]
    window = min(args.window_mb * 1_000_000, chrom_size // 2)
    window_start, _score = _choose_window(bw, chrom, chrom_size, window)
    window_end = min(chrom_size, window_start + window)
    print(
        f"chosen {chrom} (size {chrom_size}), window "
        f"{window_start}-{window_end} ({(window_end - window_start) / 1e6:.1f} Mb)"
    )

    entries = _extract_entries(
        bw, chrom, window_start, window_end, args.bin_bp, args.min_value
    )
    bw.close()
    all_bins = (window_end - window_start) // args.bin_bp
    print(
        f"entries: {len(entries)} / {all_bins} bins "
        f"(threshold {args.min_value})"
    )

    bundle_name = "ENCODE_demo_signal.bw"
    bundle_path = os.path.join(args.out_dir, bundle_name)
    summary = write_bigwig(
        bundle_path,
        [(chrom, chrom_size)],
        entries,
        compresslevel=args.compress_level,
    )
    print(
        f"bundle: {bundle_path} -> {summary['bytes']:,} bytes "
        f"({summary['n_blocks']} blocks)"
    )

    with open(bundle_path, "rb") as fh:
        bundle_sha256 = hashlib.sha256(fh.read()).hexdigest()

    meta = _experiment_metadata()
    provenance = {
        "source_file": SOURCE_ACCESSION,
        "source_url": (
            f"https://www.encodeproject.org/files/{SOURCE_ACCESSION}/"
            f"@@download/{SOURCE_ACCESSION}.bigWig"
        ),
        "source_md5": SOURCE_MD5,
        "source_size_bytes": SOURCE_SIZE_BYTES,
        "experiment": EXPERIMENT,
        "experiment_description": meta["description"],
        "target": meta["target"],
        "organism": "Mus musculus",
        "biosample": meta["biosample"],
        "assembly": ASSEMBLY,
        "output_type": "fold change over control",
        "policy": (
            "ENCODE portal public data; no usage restrictions "
            "(ENCSR000CDQ is released without embargo)"
        ),
        "subset": {
            "chromosome": chrom,
            "chrom_size": chrom_size,
            "start": window_start,
            "end": window_start + window,
            "bin_bp": args.bin_bp,
            "bin_aggregation": "mean of per-base fold change within bin",
            "min_value_threshold": args.min_value,
            "total_bins": all_bins,
            "retained_bins": len(entries),
            "retained_chromosomes": [chrom],
        },
        "subset_size_bytes": summary["bytes"],
        "bundle_sha256": bundle_sha256,
        "writer": "pychipseeker/scripts/build_bigwig_py.py (pure Python)",
    }
    json_path = os.path.join(args.out_dir, "ENCODE_demo_signal.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(provenance, fh, indent=2)
    print("provenance:", json_path)

    # ---- read-back verification via pyBigWig (the package's reader) ----
    check = pyBigWig.open(bundle_path)
    b_chroms = check.chroms()
    assert dict(b_chroms) == {chrom: chrom_size}, b_chroms
    sampled = check.intervals(chrom, window_start, window_start + 200_000) or []
    expected_here = [e for e in entries if window_start <= e[1] < window_start + 200_000]
    assert len(sampled) == len(expected_here), (
        len(sampled),
        len(expected_here),
    )
    for (s, e, v), (_c, es, ee, ev) in zip(sampled, expected_here):
        # values are stored as float32; compare with relative tolerance
        assert s == es and e == ee and abs(v - ev) <= max(1e-6, abs(ev) * 1e-5), (
            s,
            e,
            v,
            es,
            ee,
            ev,
        )
    stats_mean = check.stats(
        chrom, window_start, window_end, nBins=1, type="mean"
    )
    check.close()
    print(
        "read-back OK:",
        len(sampled),
        "first-region intervals match; window mean:",
        stats_mean,
    )
    return 0


if __name__ == "__main__":
    main()
