# ExecPlan — Figures 3–5: line-plus-dots styling (#1) + consistent method color/shape (#2)

**Repo:** `css-sims`. **Tracking issues:** `gregfaletto/cluster-stability-selection-draft#1` and `#2` (issues live in the paper repo; cross-repo refs do **not** auto-close — author closes on merge). **Branch:** `line-dots-consistent-colors`. One PR for both (they edit the same three ggplot constructors). *Revised after independent plan review — changes from v1 marked ⟲.*

---

## Background / current state (verified)

The simulation figures `fig_3_*` (Fig 3, `sim_5_1.R`), `fig_4_*` (Fig 4, `sim_5_3.R`), `fig_7_*` (Fig "5"/7, `sim_5_extra.R`) and the `sim_{1,2,3}_*_supp` variants are built by three functions in **`Helper Functions/toy_ex_slide_funcs.R`** — *not* the `…Plot4`/`…Plot3` variants in `simFunctions.R` (those have zero figure-driver callers):

| Function | Line | Panel | x, y |
|---|---|---|---|
| `createLossesPlot3` | 4692 | left | ModelSize, MSE |
| `createNSBStabPlot2` | 4739 | middle | ModelSize, NSBStability |
| `createStabMSEPlot2` | 4786 | right | NSBStability, MSE |

Shared shape:
```r
plot <- ggplot(df_gg, aes(x=…, y=…, color=Method, shape=Method)) +
    scale_shape_manual(values=1:n_methods) + …
if(line){ plot <- plot + suppressWarnings(geom_line()) }     # lines ONLY
else    { plot <- plot + suppressWarnings(geom_point(size=2.5, alpha=1)) }  # dots ONLY
```

1. **#1 (line-plus-dots):** `line` is an *either/or* switch. Current driver usage: `sim_5_1.R`/`sim_5_3.R` pass `line=TRUE` on left+middle main panels (→ **lines only**), default `line=FALSE` on right + all supplement panels (→ **dots only**). ⟲ **`sim_5_extra.R` (Fig 7) passes NO `line=` anywhere — all its panels are dots-only.** "Line-plus-dots" needs **both** geoms drawn together.
2. **#2 (consistent color):** there is **no `scale_color_manual`** — only `scale_shape_manual(values=1:n_methods)`. Both color (default hue scale) **and** shape (positional `1:n_methods`) are assigned over whatever method subset a figure plots, so the *same* method gets a different color **and** a different shape per figure.
3. ⟲ **`Method` is a character vector, not a factor** (`genPlotDf` builds it via `nameMap(...)` into a `data.frame` with `stringsAsFactors=FALSE`, `:5028`); ggplot coerces it to a factor with **alphabetical** levels over the present subset at plot time. Named-value scales match by **name**, so this is fine — but we must not rely on factor level order.

`nameMap()` (`:118`) maps system names → these **14** canonical display names (exact set, verified): `Lasso`, `Stability Selection`, `SS (Elastic Net)`, `Sparse CSS`, `Sparse CSS (est. clusts)`, `Weighted Averaged CSS`, `Weighted Averaged CSS (est. clusts)`, `Simple Averaged CSS`, `Simple Averaged CSS (est. clusts)`, `Cluster Rep. Lasso`, `Cluster Rep. Lasso (est. clusts)`, `Protolasso`, `Protolasso (est. clusts)`, `Elastic Net`.

---

## Goal

Every simulation figure renders each method as **connected points (line + dots)**, and **each method keeps one fixed color AND one fixed shape across all figures** — Stability Selection above all — regardless of which subset a figure shows.

---

## Changes

### A. `Helper Functions/toy_ex_slide_funcs.R`

**A1 — line-plus-dots (#1).** In all three functions, draw points **always**, add the connecting line when `line=TRUE`:
```r
plot <- plot + suppressWarnings(geom_point(size=2.5, alpha=1))
if(line){ plot <- plot + <connecting line> }
```
- `createLossesPlot3`, `createNSBStabPlot2` (x = ModelSize): `<connecting line>` = `suppressWarnings(geom_line())` — connects in x = ModelSize order, correct.
- `createStabMSEPlot2` (x = NSBStability): use `suppressWarnings(geom_path())` and sort `df_gg` by `Method, ModelSize` at the top, so each method's path traces in **model-size** order (a plain `geom_line()` would connect in stability order → wrong zig-zag). ⟲ Review confirmed `df_gg` here **always has `ModelSize`** (`genPlotDf:5028`); keep only a defensive `if("ModelSize" %in% names(df_gg))` guard (path if present, else dots-only) rather than a real fallback.

`line=TRUE` thus means "line + dots"; the only `line=TRUE` callers are main-text left/mid panels we *want* line+dots → no regression. `line=FALSE` stays dots-only.

**A2 — consistent color + shape (#2).** Add two shared maps near `nameMap`:
```r
methodColorMap <- function() c(
  "Lasso"="grey40", "Stability Selection"="#E41A1C", "SS (Elastic Net)"="#FB9A99",
  "Sparse CSS"="#377EB8", "Sparse CSS (est. clusts)"="#80B1D3",
  "Weighted Averaged CSS"="#4DAF4A", "Weighted Averaged CSS (est. clusts)"="#B3DE69",
  "Simple Averaged CSS"="#984EA3", "Simple Averaged CSS (est. clusts)"="#BC80BD",
  "Cluster Rep. Lasso"="#FF7F00", "Cluster Rep. Lasso (est. clusts)"="#FDB462",
  "Protolasso"="#A65628", "Protolasso (est. clusts)"="#D9A679", "Elastic Net"="#999999")
methodShapeMap <- function() c(
  "Lasso"=16,"Stability Selection"=17,"SS (Elastic Net)"=2,"Sparse CSS"=15,
  "Sparse CSS (est. clusts)"=0,"Weighted Averaged CSS"=18,"Weighted Averaged CSS (est. clusts)"=5,
  "Simple Averaged CSS"=8,"Simple Averaged CSS (est. clusts)"=7,"Cluster Rep. Lasso"=3,
  "Cluster Rep. Lasso (est. clusts)"=4,"Protolasso"=6,"Protolasso (est. clusts)"=9,"Elastic Net"=1)
```
In each function, **replace** `scale_shape_manual(values=1:n_methods)` with `scale_shape_manual(values=methodShapeMap())` and **add** `scale_color_manual(values=methodColorMap())`. Named values match by level name → each method's color+shape are fixed across every figure (subset just uses fewer keys; legend shows only present methods — no 14-way bloat). ⟲ This also removes `createNSBStabPlot2`'s reliance on the global `n_methods` for shapes.

⟲ **Coverage guard (S1):** unmapped methods otherwise render silent grey in a published figure. Add to each function, before the scales:
```r
stopifnot(all(unique(as.character(df_gg$Method)) %in% names(methodColorMap())))
```
so a missing/mistyped key fails loudly. (All figure methods come from `nameMap`'s 14, so this never fires in legitimate use.)

⟲ **Scope note (#2 → also shape).** The issue asks for consistent *color*. We also fix *shape* because fixing color alone leaves shape drifting positionally → the shared color/shape legend desyncs across figures (a new inconsistency). Doing both is the same mechanism and fully delivers "consistent method identity." Easily reverted to color-only if the author prefers; flagged in the PR.

### B. Drivers — `line=TRUE` on the main-text Figure 3/4/7 panels only

⟲ Scoped to the literal issue ("Figures 3–5") = the three **main-text** simulation figures. Add `line=TRUE` to the main-text left/mid/right panels that lack it (10 calls):

- **`sim_5_1.R`** (Fig 3; left/mid already `line=TRUE`): right panels **383 (known), 451 (est)**.
- **`sim_5_3.R`** (Fig 4; left/mid already `line=TRUE`): right panels **335 (known), 400 (est)**.
- **`sim_5_extra.R`** (Fig 7; ⟲ **nothing has `line=TRUE`**): left/mid/right of both blocks — **321, 326, 329 (known), 383, 388, 392 (est)**.

`scale_color_manual`/`scale_shape_manual` need **no** driver change — they apply to every caller automatically, so #2's color/shape consistency lands on *all* figures at once. Legend extraction differs per driver (`get_legend(fig_*_right…)` in 5_1/5_3, `fig_*_left…` in 5_extra) but is unaffected. (Line numbers pre-edit; implementer re-confirms by content.)

## Out of scope (flag in PR)

- ⟲ **Line-plus-dots is added to the main-text Figures 3/4/7.** The **intro Figure 2** MSE panel (`sim_5_1.R:341 fig_2_right`) *already* passed `line=TRUE` (lines-only), so with the function change it **also** becomes line-plus-dots automatically — a consistency win, disclosed not hidden. The **supplement twins** (`sim_{1,2,3}_*_supp`) and the **real-data figure** (`plant.R → fig_real_data`) keep their current **dots-only** styling — they aren't "Figures 3–5". All figures' **colors/shapes do** become consistent (the scale change is global, an unambiguous improvement). Making the remaining dots-only figures line-plus-dots too is a trivial follow-up (`line=TRUE` on those calls) — flagged for a one-word go-ahead.
- ⟲ **`plant.R` / `fig_real_data.pdf` is color-affected**, not untouched: it calls the same three functions, so on next run its colors switch to the fixed palette (it uses the default `line=FALSE` → no line+dots change; all its methods are mapped → no greying). Likely an improvement; **author to confirm the recolor is acceptable** at merge.
- **Regenerating the figure PDFs** in `figures/` — needs a full pipeline run (`n_sims=2000`, no cached results, `doMC`/`doParallel` etc.). This PR changes plotting *code*; figures regenerate next run (ideally batched with #6/#7/#9). Styling verified on synthetic data instead.
- ⟲ **`simFunctions.R:531` `dataResults`→`createLossesPlot3(losses_mat, …)`**: a pre-broken, figure-driver-unused caller passing a matrix (no `ModelSize/Method`); it would already error on the existing `aes`. Our scale/assert additions don't worsen it — noted so a future reader doesn't blame this change.

## Risks
- **Palette legibility** for 14 categories — mitigated by hue families + distinct shapes (the `(est. clusts)` twin shares hue family but differs in shape); colors author-adjustable. Twins rarely co-occur in one figure (known-vs-est figures are separate).
- **Shapes 0–18 distinguishability** — chosen to avoid fill-requiring shapes (21–25); fine for ≤~10 methods/panel.

## Verification (main loop, R available) ⟲ strengthened
Source `toy_ex_slide_funcs.R`; build synthetic `df_gg` (cols `ModelSize, MSE, MSELower, MSEUpper, NSBStability, StabLower, StabUpper, Method`). Render all three functions and assert:
1. **line+dots**: output has both a line layer and a point layer when `line=TRUE`.
2. **color/shape stability**: two disjoint method subsets that **both include "Stability Selection"** → SS identical color *and* shape in both; build a side-by-side PNG.
3. ⟲ **SS-omitting subset** renders fine (no crash, other methods keep their colors).
4. ⟲ **unmapped-method case** (a bogus `Method`) → must **error** at the `stopifnot`, not silently grey.
5. ⟲ **non-monotone NSBStability vs ModelSize** in `createStabMSEPlot2` → path follows ModelSize order (visual).
6. `plot_errors` both ways for Losses/NSBStab (include `*Lower/*Upper` cols); for `createStabMSEPlot2` use `plot_errors=FALSE` (⟲ its `plot_errors=TRUE` path references an undefined global `cluster_count` at `:4814`; real right-panel calls all use `FALSE`).

## Files changed
- `Helper Functions/toy_ex_slide_funcs.R` — 3 plot functions + `methodColorMap()` + `methodShapeMap()` + coverage `stopifnot`.
- `sim_5_1.R`, `sim_5_3.R`, `sim_5_extra.R` — `line=TRUE` per the §B checklist.
- `plans/issue-1-2-line-dots-consistent-colors.md` — this plan.
- (No figure PDFs regenerated — see Out of scope.)
