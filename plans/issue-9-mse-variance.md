# ExecPlan — #9: compute & report variance of test-set squared errors (rev. after review)

**Repo:** `css-sims`. **Tracking issue:** `gregfaletto/cluster-stability-selection-draft#9` (cross-repo; author closes on merge). **Branch:** `issue-9-mse-variance`.

**Coordination:** the working tree has the author's in-progress edits to `sim_5_3.R` / `sim_5_extra.R` (unequal-noise run setup). This PR touches **only** `Helper Functions/eval_functions.R` (+ a small aggregation helper) and this plan — those drivers are **not** staged. All changes are backward-compatible (trailing default arg), so a run kicked off now is unaffected, and the new metric can be `evaluate()`d on an already-run sim without re-simulating.

## The advisor ask (verbatim, `Notes from 20231204 meeting.txt:9`)
> "calculating variance of squared errors on test set (n_test observations), divide that by n_test, take the means of those values across all n_sims simulations, report that value. (take the one that is larger according to Jensen's inequality)"

**Deliverable (primary, `A`):** per sim *s*, the test-set squared errors are `e_{s,j} = (mu_hat_{s,j} - testMu_j)^2`, `j = 1..n_test`. Compute `var_j(e_{s,j})`, divide by `n_test`, then average across sims:
`A = mean_s( var_j(e_{s,j}) / n_test )`.

**On the Jensen clause (resolved interpretation):** "take the one that is larger according to Jensen's inequality" is **not** a max over two estimators — Jensen relates `f(E[X])` to `E[f(X)]`, not two variances. Its most defensible meaning is that it **justifies using *squared* errors**: with `r = mu_hat - testMu`, Jensen gives `mean_j(r^2) ≥ (mean_j r)^2` (i.e. MSE ≥ bias²), so working with the *squared* errors yields the larger/conservative quantity. Under this reading the deliverable is exactly `A` — no `B`/`max` needed. We therefore make `A` the headline and additionally expose, for the author's comparison (clearly labeled, not the answer): `B = var_s(MSE_s)` (across-sims spread of per-sim MSEs) and `B/n_sims` (to match the existing `sd/sqrt(n)` SE convention in `df_sim_stats`/`genPlotDf`).

**One question to confirm with the author/advisor (flag in PR, do not block):** in "variance of squared errors," is `e_j` the **squared error** `(mu_hat-testMu)^2` (so we take the variance *of squared errors* → `A`, the literal & Jensen-consistent reading we implement), or the **residual** `(mu_hat-testMu)` (variance of errors, squaring entering only via Jensen)? These give different numbers; we implement the former and surface the question.

## Current architecture (verified by review, `eval_functions.R`)
`cssr_mse <- new_metric("cssr_mse", ...)` (L225) dispatches on output shape to **four** branches across **three** funcs, each returning a length-`max_model_size` vector of **mean** squared errors:
- `css_res` → `cssr_mse_metric_func(out, max_model_size)` (L415): `mu_hat_i <- getCssPreds(...)`; `mses[i] <- mean((mu_hat_i - out$testMu)^2)`. Per-point errors directly available.
- `selected_clusts_list` → `clus_lasso_metric_func(out, max_model_size)` (L545) → `get_mse` (L524).
- `selected_sets` (protolasso) **and** `lasso_selected` (lasso/EN) → `lasso_metric_func(...)` (L484) → `get_mse`.

`get_mse(x, y, mu)` (L524) fits `lm`, predicts, returns `mean((preds - mu)^2)`. Per-point squared errors live in exactly two places: `get_mse` and `cssr_mse_metric_func`. **Both paths refit on the test set and predict in-sample** (`getCssPreds` is called with `trainX=testX, trainY=testY`; `get_mse` fits on `(testX[,sel], testY)`), so the squared-error variances are comparable across all method families — confirmed by review.

**Other callers of these funcs (must stay unaffected):** `weight_mse` (L382) and `cssr_mse_plant` (L346) call them positionally → the new `stat` arg **must be trailing** with default `"mean"`.

## Changes — `Helper Functions/eval_functions.R`

**1. Thread a TRAILING `stat` parameter (default `"mean"` → byte-identical existing behavior).**
- `get_mse(x_train, y_train, mu_train, stat = "mean")`: `sq <- (preds - mu_train)^2`; return `if (stat == "mean") mean(sq) else var(sq)`. (Note: returns **raw `var(sq)`**, *not* `/n_test` — see item 3.)
- `cssr_mse_metric_func(out, max_model_size, stat = "mean")`: `sq_i <- (mu_hat_i - out$testMu)^2`; `vals[i] <- if (stat == "mean") mean(sq_i) else var(sq_i)`.
- `lasso_metric_func(selected, out, max_model_size, stat = "mean")` and `clus_lasso_metric_func(out, max_model_size, stat = "mean")`: pass `stat` through to `get_mse`.
- `cssr_mse` metric: unchanged (funcs called with the default).

**2. New metric `cssr_mse_var`** = `new_metric("cssr_mse_var", "Variance of test-set squared errors", metric = function(model, out){...})` whose body **copies all four `if/else if` dispatch branches of `cssr_mse` verbatim**, passing `stat = "var_over_n"`→ actually `stat="var"` (raw). Returns, per (method, sim, model size), `var_j(e)` (raw variance of the squared errors). `NA` exactly where `cssr_mse` is `NA` (same skip logic).

**3. Aggregation helper `mseVarReport(e_df, max_model_size, n_test, n_sims)`** (new). The raw eval frame from `as.data.frame(evals(sim))` has **no `ModelSize` column** — model size is positional (`max_model_size` rows per (Method, Draw), in order), exactly as `genPlotDf` assumes via `edf_inds_k <- max_model_size*(0:(n_sims-1)) + k`. So `mseVarReport` reconstructs `ModelSize` the same way (verify the frame's row layout in the smoke harness first), then per (Method, ModelSize) computes with `na.rm = TRUE`:
- `mean_var_sq = mean_s( cssr_mse_var )`  → **`A = mean_var_sq / n_test`** (the headline reported value),
- `B = var_s( cssr_mse )` and `B_over_nsims = B / n_sims` (labeled comparison columns).
Returns a tidy frame `Method, ModelSize, A (var_reported), mean_var_sq, B, B_over_nsims`. Kept separate from `genPlotDf`/`df_sim_stats` (which stay untouched); the eventual LaTeX table consumes this frame.

## Out of scope (flag in PR)
- **Checkbox 2 — "decide where to report it (main text vs supplement tables) and update the LaTeX."** Needs real values from a full run (the author's in-progress sims) and the author's placement decision. Deferred; this PR makes the numbers computable.
- **#10 — stability metric can't use the same variance calculation.** Not addressed; the NSB-stability metric is not a test-set mean, so its uncertainty needs its own treatment (separate issue).
- No figures/tables regenerated; `genPlotDf`/`df_sim_stats` untouched.

## Risks
- **Jensen / squared-vs-residual estimand** — the one judgment call; implemented as the literal "variance of squared errors" (`A`), flagged for author confirmation (above).
- `var()` of length-1 → `NA`; only if `n_test == 1` (it is 1000–10000) — non-issue.
- `var_s` over <2 non-NA sims → `NA`; use `na.rm=TRUE`, tolerate sparse large-model-size cells.
- `stat` must be **trailing** so `weight_mse`/`cssr_mse_plant` positional calls are unaffected.

## Verification (smoke harness, NO full run)
Extend the working smoke harness (small `n_sims`, serial/few cores) to `evaluate(list(cssr_mse, cssr_mse_var))`, then:
0. **Inspect `as.data.frame(evals(sim))`** column/row layout to confirm the positional model-size stride before writing `mseVarReport`.
1. **Regression:** `cssr_mse` values are identical to a run without the change (the `stat="mean"` default path is untouched).
2. `cssr_mse_var` is finite and ≥ 0 for every (method, model size) that has a non-NA `cssr_mse`, `NA` exactly where `cssr_mse` is `NA`.
3. `mseVarReport` returns `A, B, B_over_nsims` with sane magnitudes (`A` ≈ on the order of `(SE of MSE)^2`, i.e. small; `B` on the order of the across-sims MSE spread).
4. All **four** dispatch branches exercised: SS (`css_res`), simple-avg CSS (`css_res`), cluster-rep lasso (`selected_clusts_list`), protolasso (`selected_sets`), lasso (`lasso_selected`).

## Files changed
- `Helper Functions/eval_functions.R` — trailing `stat` on `get_mse` + 3 funcs; new `cssr_mse_var` metric (all four branches); `mseVarReport` helper.
- `plans/issue-9-mse-variance.md` — this plan.
- (Drivers, figures, tables untouched; the author's `sim_5_3`/`sim_5_extra` edits not staged.)
