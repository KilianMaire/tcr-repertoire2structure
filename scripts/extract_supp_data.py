"""Extract committed CSVs for supplementary figures S3 and S4 from the raw folds.

S3: mean chain_pair_iptm (5x5) over the 5 Protenix samples of one confident
    cognate complex -> paper/data/chain_pair_iptm_example.csv
S4: per-sample iptm_TCRpep_max (= max(cpi[TCRa][pep], cpi[TCRb][pep])) for a
    cognate and its composition-scramble -> paper/data/per_sample_readouts.csv

Chain index order: 0=TCRa (A), 1=TCRb (B), 2=MHC (C), 3=b2m (D), 4=peptide (E).
Run once locally (needs runs/, GPU-produced); the CSVs it writes are committed and
every supplementary figure reads only those.
"""
from __future__ import annotations
import csv, glob, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "paper/data"
CHAINS = ["TCRa", "TCRb", "MHC", "b2m", "peptide"]

COGNATE = ROOT / ("runs/hla_a1101/folds/01685694bec8__AVFDRKSDAK/"
                  "01685694bec8__AVFDRKSDAK/seed_101/predictions")
SCRAMBLE = ROOT / ("runs/hla_a1101/folds/01685694bec8____scramble__/"
                   "01685694bec8____scramble__/seed_101/predictions")


def _samples(pred_dir):
    out = []
    for f in sorted(glob.glob(str(pred_dir / "*summary_confidence_sample_*.json"))):
        out.append(json.load(open(f)))
    return out


def _iptm_tcrpep_max(d):
    cpi = d["chain_pair_iptm"]
    return max(cpi[0][4], cpi[1][4])


def write_chain_pair_matrix():
    samples = _samples(COGNATE)
    n = len(samples)
    # mean over samples
    mat = [[sum(s["chain_pair_iptm"][i][j] for s in samples) / n for j in range(5)] for i in range(5)]
    with (DATA / "chain_pair_iptm_example.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row_chain", "col_chain", "iptm"])
        for i in range(5):
            for j in range(5):
                w.writerow([CHAINS[i], CHAINS[j], f"{mat[i][j]:.4f}"])
    print(f"chain_pair_iptm_example.csv written (mean over {n} samples)")


def write_per_sample():
    rows = []
    for construct, pred in (("cognate", COGNATE), ("scramble", SCRAMBLE)):
        for k, d in enumerate(_samples(pred)):
            rows.append({"construct": construct, "sample": k,
                         "iptm_TCRpep_max": f"{_iptm_tcrpep_max(d):.4f}",
                         "iptm_groove": f"{d['chain_pair_iptm'][2][4]:.4f}",
                         "ranking_score": f"{d['ranking_score']:.4f}"})
    with (DATA / "per_sample_readouts.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["construct", "sample", "iptm_TCRpep_max",
                                          "iptm_groove", "ranking_score"])
        w.writeheader()
        w.writerows(rows)
    print(f"per_sample_readouts.csv written ({len(rows)} rows)")


if __name__ == "__main__":
    write_chain_pair_matrix()
    write_per_sample()
