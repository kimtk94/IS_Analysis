#!/usr/bin/env python3
"""Install-free, checkpointed harmonization, MR, sensitivity and colocalization.

The input contract deliberately separates selected instruments from complete
regional summary statistics.  Colocalization never reuses the instrument-only
table.  Every stage writes a status document, including scientifically
meaningful NOT_RUN/FAILED outcomes, and may therefore be resumed independently.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (root / path).resolve()


def checkpoint(out: Path, stage: str, status: str, reason: str, **details) -> None:
    target = out / stage / "status.json"; target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({"stage": stage, "status": status, "reason": reason, **details}, indent=2) + "\n", encoding="utf-8")


def allele_state(ea: str, oa: str, target_ea: str, target_oa: str) -> tuple[str, int]:
    pair = (ea.upper(), oa.upper()); target = (target_ea.upper(), target_oa.upper())
    pal = set(pair) in ({"A", "T"}, {"C", "G"})
    if pair == target: state, sign = "DIRECT", 1
    elif pair == target[::-1]: state, sign = "SWAPPED", -1
    elif tuple(x.translate(COMPLEMENT) for x in pair) == target: state, sign = "COMPLEMENT", 1
    elif tuple(x.translate(COMPLEMENT) for x in pair) == target[::-1]: state, sign = "SWAPPED_COMPLEMENT", -1
    else: return "MISMATCH", 0
    return ("PALINDROMIC_" + state if pal else state), sign


def harmonize(cfg: dict, root: Path, out: Path) -> None:
    exposure = read_tsv(resolve(root, cfg["instrument_exposure"])); outcome = read_tsv(resolve(root, cfg["instrument_outcome"]))
    by_id = {r["variant_id"]: r for r in outcome}; rows = []
    low, high = cfg["palindromic_eaf_range"]
    for e in exposure:
        o = by_id.get(e["variant_id"])
        base = {"variant_id": e["variant_id"], "exposure_build": e["build"], "outcome_build": o["build"] if o else "",
                "exposure_chromosome": e["chromosome"], "outcome_chromosome": o["chromosome"] if o else "",
                "exposure_position": e["position"], "outcome_position": o["position"] if o else "", "allele_state": "MISMATCH",
                "eaf_difference": "", "decision": "EXCLUDE", "reason": "VARIANT_NOT_IN_OUTCOME"}
        if not o: rows.append(base); continue
        coordinates = e["build"] == o["build"] and e["chromosome"] == o["chromosome"] and e["position"] == o["position"]
        state, sign = allele_state(o["effect_allele"], o["other_allele"], e["effect_allele"], e["other_allele"])
        base["allele_state"] = state
        if not coordinates: base["reason"] = "BUILD_OR_COORDINATE_MISMATCH"
        elif state == "MISMATCH": base["reason"] = "ALLELE_MISMATCH"
        elif state.startswith("PALINDROMIC"):
            ee, oe = float(e["eaf"]), float(o["eaf"]); diff = min(abs(ee - oe), abs(ee - (1 - oe)))
            base["eaf_difference"] = f"{diff:.8g}"
            if low <= ee <= high or diff > cfg["palindromic_max_eaf_difference"]:
                base["reason"] = "PALINDROMIC_AMBIGUOUS"; rows.append(base); continue
            # Frequency evidence may reject ambiguity, but never chooses a direction:
            # only an already unambiguous allele ordering supplies the sign.
            base["reason"] = "PALINDROMIC_EXCLUDED_BY_POLICY"; rows.append(base); continue
        else:
            base.update({"decision": "INCLUDE", "reason": "ALIGNED", "beta_exposure": e["beta"], "se_exposure": e["se"],
                         "eaf_exposure": e["eaf"], "beta_outcome": f"{float(o['beta']) * sign:.12g}", "se_outcome": o["se"],
                         "eaf_outcome": o["eaf"]})
        rows.append(base)
    fields = ["variant_id", "exposure_build", "outcome_build", "exposure_chromosome", "outcome_chromosome",
              "exposure_position", "outcome_position", "allele_state", "eaf_difference", "decision", "reason",
              "beta_exposure", "se_exposure", "eaf_exposure", "beta_outcome", "se_outcome", "eaf_outcome"]
    write_tsv(out / "harmonization" / "variants.tsv", fields, rows)
    included = sum(r["decision"] == "INCLUDE" for r in rows)
    checkpoint(out, "harmonization", "SUCCESS", "VARIANT_QC_RECORDED", variants=len(rows), included=included, excluded=len(rows)-included)


def ivw(rows: list[dict], reverse: bool = False) -> tuple[float, float]:
    xkey, ykey, skey = ("beta_outcome", "beta_exposure", "se_exposure") if reverse else ("beta_exposure", "beta_outcome", "se_outcome")
    terms = [(float(r[xkey]), float(r[ykey]), 1 / float(r[skey]) ** 2) for r in rows]
    den = sum(w*x*x for x, _, w in terms); beta = sum(w*x*y for x, y, w in terms) / den
    return beta, math.sqrt(1 / den)


def normal_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2))


def mr(cfg: dict, out: Path) -> None:
    source = read_tsv(out / "harmonization" / "variants.tsv"); rows = [r for r in source if r["decision"] == "INCLUDE"]
    fields = ["method", "status", "support_condition", "instrument_count", "beta", "se", "p_value"]
    results = []
    if not rows:
        write_tsv(out / "mr" / "estimates.tsv", fields, [])
        checkpoint(out, "mr", "NOT_RUN_INSUFFICIENT_INSTRUMENTS", "NO_HARMONIZED_INSTRUMENTS", instrument_count=0); return
    if len(rows) == 1:
        r = rows[0]; beta = float(r["beta_outcome"]) / float(r["beta_exposure"])
        se = float(r["se_outcome"]) / abs(float(r["beta_exposure"])); results.append({"method": "WALD_RATIO", "status": "SUCCESS", "support_condition": "n=1", "instrument_count": 1, "beta": beta, "se": se, "p_value": normal_p(beta/se)})
    else:
        beta, se = ivw(rows); results.append({"method": "IVW", "status": "SUCCESS", "support_condition": "n>=2", "instrument_count": len(rows), "beta": beta, "se": se, "p_value": normal_p(beta/se)})
    minimum = cfg["minimum_instruments_for_robust_methods"]
    ratios = sorted((float(r["beta_outcome"])/float(r["beta_exposure"]), 1/float(r["se_outcome"])**2) for r in rows)
    for method in ("WEIGHTED_MEDIAN", "MR_EGGER", "WEIGHTED_MODE"):
        item = {"method": method, "instrument_count": len(rows), "support_condition": f"n>={minimum}"}
        if len(rows) < minimum: item["status"] = "NOT_RUN_INSUFFICIENT_INSTRUMENTS"
        else:
            if method == "WEIGHTED_MEDIAN":
                half=sum(w for _,w in ratios)/2; cumulative=0; estimate=ratios[-1][0]
                for value,w in ratios:
                    cumulative += w
                    if cumulative >= half: estimate=value; break
            elif method == "WEIGHTED_MODE": estimate=max(ratios, key=lambda x: x[1])[0]
            else:
                xs=[float(r["beta_exposure"]) for r in rows]; ys=[float(r["beta_outcome"]) for r in rows]; ws=[1/float(r["se_outcome"])**2 for r in rows]
                xm=sum(w*x for w,x in zip(ws,xs))/sum(ws); ym=sum(w*y for w,y in zip(ws,ys))/sum(ws)
                estimate=sum(w*(x-xm)*(y-ym) for w,x,y in zip(ws,xs,ys))/sum(w*(x-xm)**2 for w,x in zip(ws,xs))
            item.update(status="SUCCESS", beta=estimate, se="", p_value="")
        results.append(item)
    write_tsv(out / "mr" / "estimates.tsv", fields, results); checkpoint(out, "mr", "SUCCESS", "PRIMARY_ESTIMATE_AVAILABLE", instrument_count=len(rows))


def sensitivity(cfg: dict, out: Path) -> None:
    rows = [r for r in read_tsv(out / "harmonization" / "variants.tsv") if r["decision"] == "INCLUDE"]
    fields=["analysis", "status", "support_condition", "instrument_count", "estimate", "statistic", "p_value", "detail"]
    result=[]; n=len(rows)
    if n >= 2:
        beta,_=ivw(rows); q=sum((float(r["beta_outcome"])-beta*float(r["beta_exposure"]))**2/float(r["se_outcome"])**2 for r in rows)
        result.append({"analysis":"COCHRAN_Q","status":"SUCCESS","support_condition":"n>=2","instrument_count":n,"statistic":q,"detail":f"df={n-1}"})
    else: result.append({"analysis":"COCHRAN_Q","status":"NOT_RUN_INSUFFICIENT_INSTRUMENTS","support_condition":"n>=2","instrument_count":n})
    minimum=cfg["minimum_instruments_for_robust_methods"]
    if n >= minimum:
        xs=[float(r["beta_exposure"]) for r in rows]; ys=[float(r["beta_outcome"]) for r in rows]; ws=[1/float(r["se_outcome"])**2 for r in rows]
        xm=sum(w*x for w,x in zip(ws,xs))/sum(ws); ym=sum(w*y for w,y in zip(ws,ys))/sum(ws); slope=sum(w*(x-xm)*(y-ym) for w,x,y in zip(ws,xs,ys))/sum(w*(x-xm)**2 for w,x in zip(ws,xs)); intercept=ym-slope*xm
        result.append({"analysis":"EGGER_INTERCEPT","status":"SUCCESS","support_condition":f"n>={minimum}","instrument_count":n,"estimate":intercept})
    else: result.append({"analysis":"EGGER_INTERCEPT","status":"NOT_RUN_INSUFFICIENT_INSTRUMENTS","support_condition":f"n>={minimum}","instrument_count":n})
    if n >= 3:
        for omitted in rows:
            beta,se=ivw([r for r in rows if r is not omitted]); result.append({"analysis":"LEAVE_ONE_OUT","status":"SUCCESS","support_condition":"n>=3","instrument_count":n,"estimate":beta,"detail":"omitted="+omitted["variant_id"]})
    else: result.append({"analysis":"LEAVE_ONE_OUT","status":"NOT_RUN_INSUFFICIENT_INSTRUMENTS","support_condition":"n>=3","instrument_count":n})
    if n:
        rx=sum(2*float(r["eaf_exposure"])*(1-float(r["eaf_exposure"]))*float(r["beta_exposure"])**2 for r in rows)
        ry=sum(2*float(r["eaf_outcome"])*(1-float(r["eaf_outcome"]))*float(r["beta_outcome"])**2 for r in rows)
        result.append({"analysis":"STEIGER_DIRECTIONALITY","status":"SUCCESS","support_condition":"EAF and beta available","instrument_count":n,"estimate":rx-ry,"detail":"EXPOSURE_TO_OUTCOME" if rx>ry else "DIRECTION_NOT_SUPPORTED"})
        rb,rs=ivw(rows, reverse=True); result.append({"analysis":"REVERSE_MR","status":"SUCCESS","support_condition":"harmonized n>=1","instrument_count":n,"estimate":rb,"statistic":rs,"p_value":normal_p(rb/rs),"detail":"WALD_OR_IVW_REVERSE"})
    else:
        for name in ("STEIGER_DIRECTIONALITY","REVERSE_MR"): result.append({"analysis":name,"status":"NOT_RUN_INSUFFICIENT_INSTRUMENTS","support_condition":"harmonized n>=1","instrument_count":0})
    write_tsv(out/"sensitivity"/"results.tsv",fields,result); checkpoint(out,"sensitivity","SUCCESS" if n else "NOT_RUN_INSUFFICIENT_INSTRUMENTS","CONDITIONAL_ANALYSES_RECORDED",instrument_count=n)


def log_abf(beta: float, se: float, prior: float) -> float:
    v=se*se; return .5*(math.log(v/(v+prior)) + beta*beta*prior/(v*(v+prior)))


def coloc(cfg: dict, root: Path, out: Path) -> None:
    e=read_tsv(resolve(root,cfg["regional_exposure"])); o=read_tsv(resolve(root,cfg["regional_outcome"])); eo={r["variant_id"]:r for r in e}; oo={r["variant_id"]:r for r in o}; shared=sorted(set(eo)&set(oo)); coverage=len(shared)/max(len(e),len(o),1)
    fields=["prior_set","common_snps","exposure_region_snps","outcome_region_snps","region_coverage","pp0","pp1","pp2","pp3","pp4","status"]
    results=[]
    if len(shared) < cfg["minimum_common_snps"] or coverage < cfg["minimum_region_coverage"]:
        write_tsv(out/"colocalization"/"results.tsv",fields,[]); checkpoint(out,"colocalization","FAILED_LOW_LOCUS_COVERAGE","COMMON_SNP_OR_COVERAGE_GATE_FAILED",common_snps=len(shared),region_coverage=coverage); return
    for priors in cfg["coloc_priors"]:
        a=[math.exp(log_abf(float(eo[v]["beta"]),float(eo[v]["se"]),priors["effect_variance"])) for v in shared]; b=[math.exp(log_abf(float(oo[v]["beta"]),float(oo[v]["se"]),priors["effect_variance"])) for v in shared]
        s1,s2=sum(a),sum(b); same=sum(x*y for x,y in zip(a,b)); weights=[1,priors["p1"]*s1,priors["p2"]*s2,priors["p1"]*priors["p2"]*max(s1*s2-same,0),priors["p12"]*same]; total=sum(weights)
        results.append({"prior_set":priors["name"],"common_snps":len(shared),"exposure_region_snps":len(e),"outcome_region_snps":len(o),"region_coverage":coverage,**{f"pp{i}":weights[i]/total for i in range(5)},"status":"SUCCESS"})
    write_tsv(out/"colocalization"/"results.tsv",fields,results)
    mode=cfg.get("multiple_signal_method")
    multi_status="READY_CONDITIONAL" if mode=="conditional" else "READY_SUSIE" if mode=="susie" else "NOT_RUN_SINGLE_SIGNAL_ASSUMPTION"
    (out/"colocalization"/"multiple_signal_status.json").write_text(json.dumps({"status":multi_status,"method":mode or "single_signal"},indent=2)+"\n",encoding="utf-8")
    checkpoint(out,"colocalization","SUCCESS","FULL_REGIONAL_STATISTICS_ANALYZED",common_snps=len(shared),region_coverage=coverage,prior_sets=len(results))


def run(config_path: Path, stage: str) -> None:
    cfg=json.loads(config_path.read_text(encoding="utf-8")); root=config_path.parent; out=resolve(root,cfg["output_dir"]); out.mkdir(parents=True,exist_ok=True)
    actions={"harmonization":lambda:harmonize(cfg,root,out),"mr":lambda:mr(cfg,out),"sensitivity":lambda:sensitivity(cfg,out),"colocalization":lambda:coloc(cfg,root,out)}
    order=list(actions) if stage=="all" else [stage]
    for name in order:
        if name in ("mr","sensitivity") and not (out/"harmonization"/"variants.tsv").exists(): raise SystemExit(f"{name} requires the harmonization checkpoint")
        actions[name]()


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--config",required=True,type=Path); parser.add_argument("--stage",choices=["all","harmonization","mr","sensitivity","colocalization"],default="all"); args=parser.parse_args(); run(args.config,args.stage)
