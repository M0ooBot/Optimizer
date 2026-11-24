# Vertical fin with 5 mm cap reach

[Whole-fin STEP](vertical-fin.step) · [Outline and zebra views](fin-overview.png) · [Wrap sections](fin-sections.png) · [STEP measurements](readback.json) · [Original 40 mm cap](../README.md)

The cap reaches from Z = 240 mm to Z = 245 mm. The STEP reimports as one closed, strictly valid solid. It retains the 160 mm root chord, 85 mm chord at the cut, 65 mm leading-edge sweep at the cut, and the original 12% elliptical parent skin. The cap uses the same cubic section contraction as the 40 mm design, scaled to a 5 mm rise.

| Readback check | 40 mm cap | 5 mm cap |
| --- | ---: | ---: |
| Maximum sampled attachment gap | 0 mm | 0 mm |
| Maximum sampled attachment tangent-plane mismatch | 7.66 × 10⁻¹⁰° | 3.06 × 10⁻⁹° |
| Maximum sampled attachment normal-derivative mismatch | 2.66 × 10⁻¹¹ mm⁻¹ | 5.96 × 10⁻¹⁰ mm⁻¹ |
| Maximum attachment zebra phase-gradient mismatch | 1.59 × 10⁻¹⁰ cycles/mm | 3.56 × 10⁻⁹ cycles/mm |
| Maximum attachment zebra turn at a stripe crossing | 1.12 × 10⁻⁷° | 3.10 × 10⁻⁶° |
| Lip radius over ten measured wraps | 0.262–0.702 mm | 2.446–7.433 mm |
| Smallest radius in those sampled wrap traces | 0.262 mm | 1.88 mm |
| Near-crown principal radii | 0.710 and 49.28 mm | 7.52 and 522.14 mm |

The attachment still matches curvature to the precision of the STEP checks. The 5 mm cap bends more quickly just beyond the attachment, creating a smooth curvature peak 2.17–3.06 mm above the cut before opening out to the lip in the ten measured sections. The tightest sampled wrap radius is about 1.88 mm at 95% chord and Z = 242.17 mm. These are local radii in the sampled section planes. A process-specific radius limit was not supplied.

Both attachments and both perimeter branches were measured at 501 locations in 168 camera and stripe-angle combinations. The worst attachment zebra turn is near X = 139.02 mm, Y = −3.42 mm, Z = 240 mm. The two cap patches also meet along the shared perimeter with zero sampled gap, a maximum tangent-plane mismatch of 2.66 × 10⁻¹⁰°, and a maximum phase-gradient mismatch of 3.86 × 10⁻¹⁰ cycles/mm. Ordinary patch derivatives at the collapsed crown parameter are unavailable. The perimeter scan reaches Z = 244.999995 mm, with separate samples closer to the crown used for the radius estimate.

From the original fin study directory, run:

```powershell
python build_variant.py --cap-reach-mm 5
```

The command isolates the STEP, numerical checks, and rendered views under a folder named for the selected reach. The original study remains in the parent folder.
