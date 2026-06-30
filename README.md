# multisynth

**Synthetic control method for multiple treated units, with staggered adoption, cross-unit event-time aggregation, and pooled permutation inference.**

`multisynth` is a Stata command that extends the synthetic control method (SCM) from a single treated unit to **multiple** treated units. It runs SCM once per treated unit, aligns the per-unit treatment effects on event (relative) time, aggregates them into an overall average treatment effect on the treated (ATT), and provides pooled inference across units.

## Built on top of `synth` and `synth2`

`multisynth` is written **on top of** two existing community-contributed Stata commands:

- **`synth`** (Abadie, Diamond, and Hainmueller) — the original implementation of the synthetic control estimator. `multisynth` calls `synth` internally to fit each treated unit, so **`synth` must be installed** (`ssc install synth`).
- **`synth2`** (Yan and Chen, 2023) — an enhanced wrapper around `synth` adding placebo tests, robustness checks, and visualization. `multisynth`'s per-unit estimation engine and its Mata helpers are adapted **directly from `synth2`** and embedded inside `multisynth.ado`; consequently `synth2` itself does **not** need to be installed for `multisynth` to run.

`multisynth` is distributed as a **separate** command so that `synth2` remains intact as a reference baseline.

## What `multisynth` adds

Relative to `synth2` (which handles a single treated unit), `multisynth` adds:

1. **Multiple treated units** — `trunit()` accepts a list; the estimator runs once per treated unit.
2. **Staggered adoption** — each treated unit may have its own treatment time; effects are aligned on event time (periods relative to each unit's own treatment).
3. **Per-unit specification** — predictors and the matching/evaluation windows (`xperiod`, `preperiod`, `postperiod`, `mspeperiod`) are keyed per treated unit.
4. **Cross-unit aggregation** — per-unit effects are averaged with equal weight on event time over a common balanced window (the intersection of the per-unit event spans), producing an overall ATT.
5. **Pooled permutation inference** — a Cavallo–Galiani–Noy style permutation test across treated units, reporting two-sided, right-sided, and left-sided p-values for each event time and for the overall ATT, with placebo-donor pruning by a pre-treatment-fit cutoff.
6. **Organized output** — a compact per-unit summary, an aggregate event-time ATT table, and saved result datasets and graphs.

## Requirements

- Stata 16 or later.
- `synth` (Abadie, Diamond, and Hainmueller), available from SSC: `ssc install synth`.

## License

`multisynth` is released under the GNU General Public License v3.0 (GPL-3), consistent with the license of `synth2`, from which its estimation engine derives. See the `LICENSE` file.

## Citation

If you use `multisynth`, please cite the underlying methods and software:

- Abadie, A., and Gardeazabal, J. (2003). The Economic Costs of Conflict: A Case Study of the Basque Country. *American Economic Review*, 93(1), 113–132.
- Abadie, A., Diamond, A., and Hainmueller, J. (2010). Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California's Tobacco Control Program. *Journal of the American Statistical Association*, 105(490), 493–505.
- Abadie, A., Diamond, A., and Hainmueller, J. (2015). Comparative Politics and the Synthetic Control Method. *American Journal of Political Science*, 59(2), 495–510.
- Cavallo, E., Galiani, S., Noy, I., and Pantano, J. (2013). Catastrophic Natural Disasters and Economic Growth. *Review of Economics and Statistics*, 95(5), 1549–1561.
- Yan, G., and Chen, Q. (2023). synth2: Synthetic control method with placebo tests, robustness test, and visualization. *The Stata Journal*, 23(3), 597–624.

## Author

Ivan Zapata — Texas Tech University — izapatad@ttu.edu
