#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
External validation of pLIN v4.1 on independent hospital datasets, as registered in
output/backbone_v41/PREREGISTRATION_external_v4.1.md (branch preregistration-v4, commit dda3947).

  E1  PRJNA981541  multispecies NDM-5 hospital outbreak, USA, 2021 to 2023
  E2  PRJNA924056  endemic IMP-4 dissemination, Melbourne (Macesic et al. 2023); no public plasmid, not analysed
  E3  PRJNA1010831 carbapenemase-producing Enterobacterales, 30 UK hospital laboratories, 2021 to 2023
      (replacement for E2, Deviation 2)

Stages (each cached): download plasmid records -> exclude development accessions -> pLIN v4.1 (one session
per dataset, identical database plasmids treated as absent) -> MOB-suite, pling -> AMRFinderPlus
(carbapenemase carriers) -> alignment truth (registered design) -> endpoints.

Usage:
  python v41_external.py [--index-dir DIR]
Output: output/backbone_v41/external/
"""

import argparse
import glob
import json
import os
import random
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_kmers import adaptive_sketch, min_containment
from plin_v4 import seq_hash
from plin_v41_typer import PlinV41Release
from v4_evaluate import N_BOOT, SEED, metrics
from v4_truth_pairs import jaccard_matrix, sample_pairs
from validate_alignment_backbone import align_pair

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
OUT = os.path.join(V41, "external")
DATASETS = {"E1": ("PRJNA981541", "NDM-5"), "E2": ("PRJNA924056", "IMP-4"),
            "E3": ("PRJNA1010831", None)}                 # E3 replaces E2 (Deviation 2); None = any carbapenemase allele
CARBAPENEMASE = ("blaKPC", "blaNDM", "blaOXA-48", "blaOXA-181", "blaOXA-232", "blaOXA-244", "blaIMP", "blaVIM")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
MOB = os.path.expanduser("~/miniconda3/envs/mob_suite_env/bin/mob_typer")
PLING = os.path.expanduser("~/miniconda3/envs/pling_env/bin/pling")
AMRF = os.path.expanduser("~/miniforge3/bin/amrfinder")
MIN_LEN, MIN_N = 1000, 20


def eutil(name, **params):
    url = EUTILS + name + ".fcgi"
    body = urllib.parse.urlencode(params).encode()
    for k in range(5):
        try:
            data = urllib.request.urlopen(urllib.request.Request(url, data=body), timeout=180).read()
            time.sleep(0.4)
            return data
        except Exception:
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"E-utilities failed: {url[:120]}")


def bare(acc):
    a = str(acc).split(".")[0]
    return a[3:] if a.startswith("NZ_") else a


def download(ds, project):
    """Plasmid records linked to the BioProject (directly or through its assemblies)."""
    d = os.path.join(OUT, ds)
    meta_f = os.path.join(d, "records.tsv")
    if os.path.exists(meta_f) and os.path.getsize(meta_f) > 1:
        return pd.read_csv(meta_f, sep="\t")
    os.makedirs(os.path.join(d, "fasta"), exist_ok=True)
    ids = set(json.loads(eutil("esearch", db="nuccore", term=project, retmax=10000, retmode="json"))["esearchresult"]["idlist"])
    asm = json.loads(eutil("esearch", db="assembly", term=project, retmax=10000, retmode="json"))["esearchresult"]["idlist"]
    for k in range(0, len(asm), 100):
        link = json.loads(eutil("elink", dbfrom="assembly", db="nuccore", id=",".join(asm[k:k + 100]), retmode="json"))
        for ls in link.get("linksets", []):
            for db in ls.get("linksetdbs", []):
                ids.update(db.get("links", []))
    ids = sorted(ids)
    rows = []
    for k in range(0, len(ids), 200):
        summ = json.loads(eutil("esummary", db="nuccore", id=",".join(ids[k:k + 200]), retmode="json"))["result"]
        for uid in summ.get("uids", []):
            r = summ[uid]
            rows.append({"uid": uid, "accession": r.get("accessionversion"), "title": r.get("title", ""),
                         "length": int(r.get("slen", 0)), "biosample": r.get("biosample", ""),
                         "topology": r.get("topology", "")})
    m = pd.DataFrame(rows, columns=["uid", "accession", "title", "length", "biosample", "topology"]).drop_duplicates("accession")
    # completeness = GenBank topology "circular" (Deviation 1: the plasmid records are WGS-type records whose
    # titles never say "complete sequence")
    keep = m.title.str.contains("plasmid", case=False) & (m.topology == "circular") & (m.length >= MIN_LEN)
    m = m[keep].sort_values("accession").reset_index(drop=True)
    for k in range(0, len(m), 100):
        fa = eutil("efetch", db="nuccore", id=",".join(m.uid.astype(str).iloc[k:k + 100]), rettype="fasta", retmode="text").decode()
        for rec in fa.split(">")[1:]:
            acc = rec.split()[0]
            open(os.path.join(d, "fasta", f"{acc}.fasta"), "w").write(">" + rec)
    m["fasta"] = [os.path.join(d, "fasta", f"{a}.fasta") for a in m.accession]
    m = m[np.array([os.path.exists(f) for f in m.fasta], dtype=bool)]
    m.to_csv(meta_f, sep="\t", index=False)
    return m


def development_accessions():
    acc = set()
    tabs = [os.path.join(BASE_DIR, "output", "backbone_v4", "eval_plasmids.tsv"),
            os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv")] + \
        glob.glob(os.path.join(V41, "confirm", "*_plasmids.tsv"))
    for t in tabs:
        df = pd.read_csv(t, sep="\t", dtype=str)
        col = "plasmid_id" if "plasmid_id" in df else "accession"
        acc |= {bare(x.replace("RefSeq_", "")) for x in df[col]}
    for lst in glob.glob(os.path.join(BASE_DIR, "output", "comparator_benchmark", "*_inputs.txt")):
        acc |= {bare(os.path.basename(l.strip()).rsplit(".", 1)[0]) for l in open(lst) if l.strip()}
    return acc


def type_plin(ds, m, rel, rows_by_hash, row_of_bare):
    out = os.path.join(OUT, ds, "plin_codes.tsv")
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t", dtype=str)
    seqs = {a: "".join(str(r.seq) for r in SeqIO.parse(f, "fasta")).upper() for a, f in zip(m.accession, m.fasta)}
    fams = rel.families([(a, seqs[a]) for a in m.accession], tempfile.mkdtemp(prefix="plin_ext_"), threads=8)
    sess = rel.hx.session()
    s_sk, s_sc, s_codes, rows = [], [], [], []
    for a in m.accession:
        excl = {r for x in rows_by_hash.get(str(np.uint64(seq_hash(seqs[a]))), []) for r in row_of_bare.get(bare(x), [])}
        excl |= set(row_of_bare.get(bare(a), []))
        sk, sc = adaptive_sketch(seqs[a])
        extra = None
        if s_codes:
            extra = {"codes": np.array(s_codes), "kmin": np.array([min_containment(sk, x, sc, y) for x, y in zip(s_sk, s_sc)])}
        code, _ = sess.assign(fams[a], sk, sc, rel.n, extra=extra, add=True, exclude=excl)
        s_sk.append(sk); s_sc.append(sc); s_codes.append(code)
        rows.append({"accession": a, "pLIN_v41": ".".join(map(str, code)), "excluded_database_rows": len(excl)})
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return df.astype(str)


def run_mob(ds, m):
    d = os.path.join(OUT, ds, "mob_typer")
    os.makedirs(d, exist_ok=True)

    def one(args):
        a, f = args
        o = os.path.join(d, f"{a}.txt")
        if not os.path.exists(o):
            subprocess.run([MOB, "-i", f, "-o", o, "-n", "1"], capture_output=True)
        return a
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(one, zip(m.accession, m.fasta)))
    clean = lambda x: None if pd.isna(x) or x in ("-", "") else x
    lab = {}
    for a in m.accession:
        r = pd.read_csv(os.path.join(d, f"{a}.txt"), sep="\t", dtype=str).iloc[0]
        lab[a] = (clean(r.get("primary_cluster_id")), clean(r.get("secondary_cluster_id")))
    return lab


def run_pling(ds, m):
    work = os.path.join(OUT, ds, "pling")
    os.makedirs(work, exist_ok=True)
    open(os.path.join(work, "inputs.txt"), "w").write("\n".join(m.fasta) + "\n")
    typing = glob.glob(os.path.join(work, "out", "dcj_thresh_*_graph", "objects", "typing.tsv"))
    if not typing:
        env = dict(os.environ, PATH=os.path.dirname(PLING) + os.pathsep + os.environ["PATH"])
        subprocess.run([PLING, "cluster", "align", "inputs.txt", "out", "--sourmash", "--cores", "8", "--visualisation", "none"],
                       cwd=work, env=env, capture_output=True)
        typing = glob.glob(os.path.join(work, "out", "dcj_thresh_*_graph", "objects", "typing.tsv"))
    if not typing:
        return None
    sub = dict(pd.read_csv(typing[0], sep="\t", dtype=str).iloc[:, :2].values)
    com = dict(pd.read_csv(os.path.join(work, "out", "containment", "containment_communities", "objects", "communities.tsv"),
                           sep="\t", dtype=str).iloc[:, :2].values)
    return sub, com


def run_amr(ds, m):
    d = os.path.join(OUT, ds, "amrfinder")
    os.makedirs(d, exist_ok=True)

    def one(args):
        a, f = args
        o = os.path.join(d, f"{a}.tsv")
        if not os.path.exists(o):
            subprocess.run([AMRF, "-n", f, "-o", o, "--threads", "1"], capture_output=True)
        return a
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(one, zip(m.accession, m.fasta)))
    genes = {}
    for a in m.accession:
        o = os.path.join(d, f"{a}.tsv")
        genes[a] = set(pd.read_csv(o, sep="\t").iloc[:, 5].astype(str)) if os.path.exists(o) and os.path.getsize(o) else set()
    return genes


def truth_pairs(ds, m):
    out = os.path.join(OUT, ds, "truth_pairs.tsv")
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t")
    work = os.path.join(OUT, ds, "sketches")
    os.makedirs(work, exist_ok=True)
    ids, mat = jaccard_matrix(sorted(m.fasta), ds, work)
    rows = sample_pairs(ids, mat, 600, random.Random(SEED))
    fasta = dict(zip(m.accession, m.fasta))
    length = {a: len(next(SeqIO.parse(fasta[a], "fasta")).seq) for a in m.accession}
    tasks, order = [], []
    for r in rows:
        a, b = (r["p1"], r["p2"]) if length[r["p1"]] <= length[r["p2"]] else (r["p2"], r["p1"])
        order.append((a, b, r))
        tasks.append((fasta[a], fasta[b], length[a], length[b], [], [], False))
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(align_pair, tasks, chunksize=4))
    out_rows = []
    for (a, b, r), x in zip(order, res):
        out_rows.append({"dataset": ds, "plasmid_shorter": a, "plasmid_longer": b, "jaccard": r["jaccard"],
                         "jaccard_bin": r["jaccard_bin"], "weight": r["weight"], "AF_shorter": x["a"]["AF"],
                         "AF_min": min(x["a"]["AF"], x["b"]["AF"]), "ANI_aln": x["ANI_aln"],
                         "bb_AF_shorter": x["a"]["bb_AF"]})
    df = pd.DataFrame(out_rows)
    df["same_lineage"] = (df.AF_min >= 0.8) & (df.ANI_aln >= 99)
    df["related_backbone"] = df.bb_AF_shorter.fillna(df.AF_shorter) >= 0.5
    df.to_csv(out, sep="\t", index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index-dir", default=None)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    dev = development_accessions()
    rel = PlinV41Release(os.path.join(V41, "release"), index_dir=args.index_dir)
    h = pd.read_csv(os.path.join(rel.dir, "plasmid_hashes.tsv.gz"), sep="\t", dtype=str)
    rows_by_hash = h.groupby("hash").accession.apply(list).to_dict()
    row_of_bare = {}
    for k, a in enumerate(rel.accessions):
        row_of_bare.setdefault(bare(a.replace("RefSeq_", "")), []).append(k)

    labels, pairs, summary, genes_all = {}, [], {}, {}
    for ds, (project, gene) in DATASETS.items():
        m = download(ds, project)
        n_all = len(m)
        m = m[~m.accession.map(bare).isin(dev)].reset_index(drop=True)
        summary[ds] = {"bioproject": project, "plasmid_records": n_all, "after_excluding_development": len(m)}
        print(ds, summary[ds], flush=True)
        if len(m) < MIN_N:
            summary[ds]["analysed"] = False
            continue
        summary[ds]["analysed"] = True
        codes = type_plin(ds, m, rel, rows_by_hash, row_of_bare).set_index("accession").pLIN_v41
        summary[ds]["in_database_identical"] = int((type_plin(ds, m, rel, rows_by_hash, row_of_bare)
                                                     .excluded_database_rows.astype(int) > 0).sum())
        mob = run_mob(ds, m)
        pl = run_pling(ds, m)
        genes_all[ds] = run_amr(ds, m)
        lab = {f"pLIN v4.1 L{k}": {a: ".".join(codes[a].split(".")[:k]) for a in m.accession} for k in range(1, 7)}
        lab["MOB-suite primary"] = {a: mob[a][0] for a in m.accession}
        lab["MOB-suite secondary"] = {a: mob[a][1] for a in m.accession}
        if pl:
            fid = {os.path.basename(f).rsplit(".", 1)[0]: a for a, f in zip(m.accession, m.fasta)}
            lab["pling subcommunity"] = {a: pl[0].get(a) for a in m.accession}
            lab["pling community"] = {a: pl[1].get(a) for a in m.accession}
        labels[ds] = lab
        pairs.append(truth_pairs(ds, m))
    pairs = pd.concat(pairs, ignore_index=True)

    # endpoints, pooled over datasets; bootstrap resamples plasmids within each dataset
    rng = np.random.default_rng(SEED)
    mult = np.ones((N_BOOT, len(pairs)))
    for ds in pairs.dataset.unique():
        sel = (pairs.dataset == ds).values
        ids = sorted(set(pairs.plasmid_shorter[sel]) | set(pairs.plasmid_longer[sel]))
        pos = {i: k for k, i in enumerate(ids)}
        counts = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=N_BOOT)
        mult[:, sel] = counts[:, [pos[x] for x in pairs.plasmid_shorter[sel]]] * counts[:, [pos[x] for x in pairs.plasmid_longer[sel]]]
    methods = sorted(set.intersection(*[set(l) for l in labels.values()]))
    w = pairs.weight.values
    rows, boots = [], {}
    for meth in methods:
        pred = np.array([labels[d][meth].get(a) is not None and labels[d][meth].get(a) == labels[d][meth].get(b)
                         for d, a, b in zip(pairs.dataset, pairs.plasmid_shorter, pairs.plasmid_longer)])
        for t in ("same_lineage", "related_backbone"):
            y = pairs[t].values
            for scope in ["pooled"] + list(pairs.dataset.unique()):
                s = np.ones(len(pairs), bool) if scope == "pooled" else (pairs.dataset == scope).values
                p, r, f = metrics(pred[s], y[s], w[s])
                rows.append({"scope": scope, "method": meth, "truth": t, "precision_w": p, "recall_w": r, "F1_w": f,
                             "positive_pairs": int(y[s].sum()), "pairs": int(s.sum())})
            W = mult * w
            tp, pp, pos_ = (W * (pred & y)).sum(1), (W * pred).sum(1), (W * y).sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                pr, rc = tp / pp, tp / pos_
                boots[(meth, t)] = np.nan_to_num(np.where(pr + rc > 0, 2 * pr * rc / (pr + rc), 0.0))
    res = pd.DataFrame(rows)
    pooled = res.scope == "pooled"
    res.loc[pooled, "F1_w_lo"] = [np.percentile(boots[(m, t)], 2.5) for m, t in zip(res.method[pooled], res.truth[pooled])]
    res.loc[pooled, "F1_w_hi"] = [np.percentile(boots[(m, t)], 97.5) for m, t in zip(res.method[pooled], res.truth[pooled])]
    res.to_csv(os.path.join(OUT, "external_metrics.tsv"), sep="\t", index=False)

    def endpoint(m1, m2, t):
        f = lambda mm: float(res[(res.scope == "pooled") & (res.method == mm) & (res.truth == t)].F1_w.iloc[0])
        d = boots[(m1, t)] - boots[(m2, t)]
        return {"F1": f(m1), "F1_comparator": f(m2), "difference": f(m1) - f(m2),
                "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5)),
                "F1_lo": float(np.percentile(boots[(m1, t)], 2.5)), "F1_hi": float(np.percentile(boots[(m1, t)], 97.5))}
    ep = {"datasets": summary, "pairs": int(len(pairs)), "same_lineage_pairs": int(pairs.same_lineage.sum()),
          "related_backbone_pairs": int(pairs.related_backbone.sum()),
          "P1_lineage_L5_vs_MOBsecondary": endpoint("pLIN v4.1 L5", "MOB-suite secondary", "same_lineage"),
          "P2_backbone_L1_vs_MOBprimary": endpoint("pLIN v4.1 L1", "MOB-suite primary", "related_backbone")}

    # secondary: carbapenemase-carrying plasmids sharing pLIN levels / comparator clusters
    sec = {}
    for ds, (project, gene) in DATASETS.items():
        if ds not in labels:
            continue
        if gene:
            carriers = [a for a, g in genes_all[ds].items() if any(x.startswith("bla" + gene) for x in g)]
            prs = [(a, b) for i, a in enumerate(carriers) for b in carriers[i + 1:]]
        else:                                             # pairs carrying the same carbapenemase allele
            allele = {a: {x for x in g if x.startswith(CARBAPENEMASE)} for a, g in genes_all[ds].items()}
            carriers = [a for a in allele if allele[a]]
            prs = [(a, b) for i, a in enumerate(carriers) for b in carriers[i + 1:] if allele[a] & allele[b]]
        sec[ds] = {"gene": gene or "same carbapenemase allele", "carrier_plasmids": len(carriers), "carrier_pairs": len(prs)}
        for meth in ("pLIN v4.1 L1", "pLIN v4.1 L5", "pLIN v4.1 L6", "MOB-suite primary", "MOB-suite secondary"):
            l = labels[ds][meth]
            sec[ds][meth] = round(100 * np.mean([l[a] is not None and l[a] == l[b] for a, b in prs]), 1) if prs else None
    ep["carbapenemase_plasmids_grouped_pct"] = sec
    json.dump(ep, open(os.path.join(OUT, "external_endpoints.json"), "w"), indent=1)
    pd.set_option("display.width", 250)
    print(res[res.scope == "pooled"].pivot_table(index="method", columns="truth", values="F1_w").round(3).to_string())
    print(json.dumps({k: v for k, v in ep.items() if k != "datasets"}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
