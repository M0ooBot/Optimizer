# Rounded tip skill retest

The designer accepted the rebuilt models in all five cases. The revised skill caught the unwanted curvature reversal on the 55° wing. The fresh build also made the edge considerably sharper. The control remains available for comparison. [Whole-model comparison viewer](../review.html) · [Swept section radius plot](swept-section-radii.png) · [Paired measurements](comparison.json)

Four fresh Sol builders started without chat history and received the revised skill, common project resources and dimensional briefs. They had no access through their supplied study directories to the previous trial designs. Parent profiles and tip outlines were locked. The twisted-wing brief included its analytic parent-section definition. One fin builder produced both reaches. The five existing STEP files served as historical controls. A separate reviewer checked the delivered STEP hashes and independently reviewed the measurements.

| Model | Findings on common physical sections |
| --- | --- |
| [55° swept wing](swept55/README.md) | Three of nine control cuts had curvature reversals. The rebuild had zero on those nine, and its tightest turn occurred at the lip in each. All nine tightest radii became smaller. At 97% chord the control turned down to 0.042106 mm before opening to a 0.171915 mm lip. The rebuilt lip itself was 0.000758 mm in that section. |
| [Straight wing](flat/README.md) | Both versions rolled into the lip without a measured reversal on ten cuts. Every minimum radius grew by about 13.8%. Midchord went from 1.266 to 1.440 mm. |
| [Twisted wing](twisted/README.md) | Both versions showed the lower-side inflection associated with the cambered parent. The rebuild rolled more tightly on every common cut. At 99% chord the smallest radius fell from 0.00513 to 0.00289 mm. Its single wrap surface includes the lip in its interior. |
| [40 mm fin](fin/40/README.md) | The two STEP exports agreed on the ten measured section radii, from 0.262 to 0.702 mm. The retained height was 280 mm. |
| [5 mm fin](fin/5/README.md) | The 245 mm fin retained its shallow outline. The smallest sampled radius over all ten cuts was 1.881 mm before and 1.889 mm after. Eight individual cuts tightened. At 45% chord the minimum changed from 3.394 to 1.889 mm. The tightest spot moved to the lip on all ten cuts, and a smaller interior curvature peak remained in central sections. |

### Paired verification

All ten exports reimported as single closed, strictly valid solids. Direct surface extrema recovered the specified spans and heights. Across all pairs the maximum difference in retained-skin ray intersections was below 0.00000025 mm. The pairwise checks include mirrored geometry, attachment separation, tangent-plane angle, the surface curvature tensor, and gradients of the same diagnostic reflection field across each side of the join. The rebuilt attachment tangent mismatch stayed below 0.0000001° in this sweep. Geometry-specific builder reports also locate actual stripe transitions using their own camera sets and are linked above.

The [paired auditor](audit.py) intersects each model with the same case-specific physical planes. It evaluates signed curvature from CAD surface derivatives along each section curve on both skins of the positive cap or the full fin. On the 55° wing the planes are parallel to the leading-edge sweep, anchored at nine cut-chord fractions from 2% to 97%. Peak curvature is sampled at 401 parameters per CAD branch and the largest sampled peak is refined locally. Bracketed sign reversals are refined with a numerical root solve. The audit also compares the physical lip point with the design coordinate on every cut. It excludes collapsed endpoints. Radii are sectional values and the cuts do not establish a whole-cap lower bound.

This is one independent build per brief and skill version. The builders chose different wrap constructions. These results show what the revised wording allowed and what it failed to prevent in this trial. No manufacturing radius or allowable reduction from the control was specified.

To reproduce from the project root, run the two versions of each case, then regenerate the summary and viewer.

```powershell
foreach ($case in 'flat','twisted','swept55','fin/40','fin/5') {
    foreach ($version in 'control','rebuilt') {
        python runs/cap-trials/retest/audit.py $case $version
    }
}
python runs/cap-trials/retest/compare.py
python runs/cap-trials/build_viewer.py
```
