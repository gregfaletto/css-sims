# ExecPlan — finish line-plus-dots styling for ALL figures (completes #1)

**Repo:** `css-sims`. **Tracking:** `gregfaletto/cluster-stability-selection-draft#1` (line-plus-dots) / `#2` (colors). **Branch:** `figure-styling-supplements`.

## Why now
PR #1 scoped `line=TRUE` to the main-text Figure 3/4/7 panels and deferred the supplement panels + `plant.R`. Two things make finishing it timely: (a) the Section 5 restructure (paper #44, merged) **moved figures into the supplement**, so those panels now matter more; (b) the data-application figure is now **main-text Fig 4, sitting next to the line-plus-dots Fig 3** — so `plant.R` must match. #2's color/shape maps are global (in the plot functions), so they already apply everywhere; this PR is purely the remaining `line=TRUE`.

## The change — add `line=TRUE` to the 24 figure calls that lack it
The connecting-line capability is already built and reviewed (PR #1); every figure uses `x=ModelSize` (or sorts by it), so the geoms are valid. This just enables line-plus-dots on the remaining calls:
- **sim_5_1.R** (6): `fig_3_supp_{left,mid,right}` at L416/423/429 (known) and L484/491/497 (est).
- **sim_5_3.R** (6): `fig_4_supp_{left,mid,right}` at L367/374/380 (known) and L431/438/444 (est).
- **sim_5_extra.R** (6): `fig_4_supp_{left,mid,right}` at L360/367/373 (known) and L423/430/436 (est).
- **plant.R** (6): main data-app `fig_4_{left,mid,right}` at L308/313/317 **and** supplement `fig_4_supp_{left,mid,right}` at L350/356/362.

Each is a `createLossesPlot3`/`createNSBStabPlot2`/`createStabMSEPlot2` call; add `, line=TRUE` as the final argument (before the closing `)`), exactly as the existing main-text calls do. Touch nothing else — no method lists, params, `n_methods`, or `evaluate()` calls (those are the author's active driver setup).

## Safety
- The three plot functions already branch on `line` (geom_point always; `+ geom_line`/`geom_path` when `line=TRUE`) — verified by render in PR #1. No function changes here.
- All driver `results_df`s carry `ModelSize` (the existing dots-only figures already use `aes(x=ModelSize)`), so `geom_line` (Losses/NSBStab) and the ModelSize-sorted `geom_path` (StabMSE) are well-defined for every call, including `plant.R`.

## Out of scope
- css-sims #2 item 1 (drop the vestigial `n_methods` across ~29 sites) — separate deferred item; not bundled here.
- Regenerating the figure PDFs — needs the full run; this is code-only so the re-run styles every figure.

## Verification (no full sims)
- `Rscript -e parse()` on all four edited drivers (syntax).
- Synthetic-data render (reuse the PR #1 harness): call each function with `line=TRUE` in the supplement/plant argument style and confirm a line/path layer + a point layer; confirm `plant.R`'s `createStabMSEPlot2(..., log_mse=TRUE, line=TRUE)` builds (geom_path under log-y).
- Re-grep: every figure call in all four drivers now has `line=TRUE` (0 remaining without).

## Files changed
- `sim_5_1.R`, `sim_5_3.R`, `sim_5_extra.R`, `plant.R` — `line=TRUE` on the 24 listed calls.
- `plans/figure-styling-supplements.md` — this plan.
