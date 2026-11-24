"""Summarize the hash-bound paired section audits and plot the swept-wing result."""
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

from audit import CASES, HERE, ROOT


def main():
    summary={}
    for case,(control,_,_) in CASES.items():
        runs=[]
        for generation,path in (("control",ROOT/control),("rebuilt",HERE/case/"wing.step")):
            report=json.loads((HERE/"measurements"/case/(generation+".json")).read_text())
            if hashlib.sha256(path.read_bytes()).hexdigest()!=report["step_sha256"]:
                raise ValueError(f"Stale audit for {path}")
            runs.append(report)
        original,revised=runs
        cuts=[]
        for a,b in zip(original["sections"],revised["sections"],strict=True):
            if a["fraction"]!=b["fraction"]:
                raise ValueError("Section locations differ")
            cuts.append(dict(chord_fraction=a["fraction"],**{
                version:dict(minimum_mm=r["smallest_radius"],lip_mm=r["lip_radius"],
                             min_xyz_mm=r["at"],branches_with_reversal=sum(
                                 bool(branch["reversals"]) for branch in r["branches"]))
                for version,r in (("control",a),("rebuilt",b))}))
        summary[case]=dict(step_sha256=[r["step_sha256"] for r in runs],
            extents_mm=[r["extent_mm"] for r in runs],
            max_delta_retained_skin_mm=float(np.max(np.abs(np.asarray(original["parent"])-revised["parent"]))),
            sections=cuts)
    (HERE/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    data=summary["swept55"]["sections"]
    chord=[100*c["chord_fraction"] for c in data]
    fig,ax=plt.subplots(figsize=(8.5,4.7),layout="constrained")
    for generation,key,label,style,color in [
        ("control","lip_mm","Control · radius at lip","o--","#819098"),
        ("control","minimum_mm","Control · tightest radius","o-","#355374"),
        ("rebuilt","minimum_mm","Rebuilt · tightest radius, at lip","s-","#ad4930"),
    ]:
        ax.semilogy(chord,[c[generation][key] for c in data],style,color=color,label=label,lw=2,ms=5)
    ax.set(xlabel="Chord position at cap attachment (%)",ylabel="Section radius (mm, log scale)",
           title="55° swept wing — same physical cuts through both STEP files",xticks=chord,
           xlim=(-1,100),ylim=(.0005,2))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda n,_:f"{n:g}"))
    ax.grid(True,which="major",alpha=.22);ax.legend(frameon=False,loc="lower left")
    fig.savefig(HERE/"swept-section-radii.png",dpi=160);plt.close(fig)


if __name__=="__main__":
    main()
