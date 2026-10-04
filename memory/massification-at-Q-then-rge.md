---
name: massification-at-q-then-rge
description: "User rule - massification must be evaluated at Q=sqrt(s_bbH) and moved to the chosen scale via exact RGE flow, validated against the exact amplitude"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 660413b7-bdbc-4081-adcf-1fd98bdd9ba7
  modified: 2026-10-02T18:17:33.966Z
---

Always evaluate the massification approximation at Q = sqrt(s_bbH) (times etascfact), then move to the chosen scale (muR) via the exact RGE flow, as done in MiNNLOPS_res/ttH_NLO/virtual.f (fromQtomuR block). Never shortcut by massifying directly at muR. Check the flow against the exact amplitude.

**Why:** user's explicit instruction (2026-10-02) after I proposed massifying at muR to avoid the running.
**How to apply:** for any massification variant (bbH, ttH, other approximations), keep scale Q for the approximation and use the massive amplitude's poles for the flow; validate with an exact-amplitude check. Related: [[m1m0-massification-status]].
