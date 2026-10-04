---
name: lambdascvar-mapping
description: User rules for the massless-mapping variation in bbHyt2 MAapprox - parameter named lambdascvar, must be IR safe (lambda>0)
metadata:
  type: feedback
---

The massive->massless mapping variation is the parameter `lambdascvar` (share of m_b^2 moved into the b pT: pT'^2 = pT^2 + lambda m^2, pz kept; gg in lab frame, q qbar in partonic frame; central 1, variations 1/2 and 2 as reweighting weights). The user called it "the new tau" but wants it named lambdascvar.

**Why:** user (2026-10-03): "you cannot use a lambda that is infrared unsafe for the massless amplitude" — lambda=0 (old MATRIX q qbar map 1) lets pT'->0, beam-collinear singular; that was also why q qbar MA undershot at low pT (0.70-0.83 below pT_b 5 GeV).
**How to apply:** never use lambda <= 0 (code exits); MAapprox_mapping 2 selects the family. Variations only on the HT/4 run unless asked. Related: [[m1m0-massification-status]], [[massification-at-Q-then-rge]].
