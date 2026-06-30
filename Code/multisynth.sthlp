{smcl}
{* *! version 0.6.0  29jun2026}{...}
{cmd:help multisynth}
{hline}
{vieweralsosee "" "--"}{...}
{vieweralsosee "synth2" "help synth2"}{...}
{vieweralsosee "synth" "help synth"}{...}
{viewerjumpto "Syntax" "multisynth##syntax"}{...}
{viewerjumpto "Description" "multisynth##description"}{...}
{viewerjumpto "Required settings" "multisynth##required"}{...}
{viewerjumpto "Options" "multisynth##options"}{...}
{viewerjumpto "Output files" "multisynth##outputs"}{...}
{viewerjumpto "Stored results" "multisynth##results"}{...}
{viewerjumpto "References" "multisynth##references"}{...}
{viewerjumpto "Author" "multisynth##author"}{...}

{title:Title}

{phang}
{bf:multisynth} {hline 2} Synthetic control method for multiple treated units, with staggered adoption and event-time aggregation

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:multisynth} {depvar}{cmd:,}
{opt tru:nit(numlist)}
{opt trp:eriod(unit: # [|| unit: # ...])}
{opt pred:ictors(unit: varlist [|| unit: varlist ...])}
[{it:options}]

{synoptset 34 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt tru:nit(numlist)}}list of treated unit numbers (one or more){p_end}
{synopt:{opt trp:eriod(unit: # || ...)}}treatment time for each treated unit (keyed){p_end}
{synopt:{opt pred:ictors(unit: varlist || ...)}}predictor list for each treated unit (keyed){p_end}

{syntab:Per-unit periods (keyed; optional)}
{synopt:{opt xp:eriod(unit: numlist || ...)}}covariate-averaging periods, per unit{p_end}
{synopt:{opt prep:eriod(unit: numlist || ...)}}pretreatment periods, per unit{p_end}
{synopt:{opt postp:eriod(unit: numlist || ...)}}posttreatment periods, per unit{p_end}
{synopt:{opt mspep:eriod(unit: numlist || ...)}}MSPE-minimization periods, per unit{p_end}

{syntab:Donor pool and output}
{synopt:{opth cou:nit(numlist)}}override the donor pool (default: never-treated units){p_end}
{synopt:{opt saver:esults(folder)}}output folder (default: {bf:multisynth_results}){p_end}
{synopt:{opt det:ail}}show full per-unit {bf:synth2} output (default: compact){p_end}
{synopt:{opt nofig:ure}}do not build any graphs{p_end}

{syntab:Passed through to the per-unit engine}
{synopt:{it:synth2 options}}e.g. {opt nested}, {opt allopt}, {opt placebo()}, {opt loo}, {opt customV()}, {opt sigf()}, ...{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}{helpb xtset} {it:panelvar} {it:timevar} must be used to declare a balanced panel in the usual long form; see {manhelp xtset XT:xtset}.{p_end}
{p 4 6 2}Only {depvar} appears before the comma; all covariates go in {opt predictors()}. Treatment times and predictors are specified {bf:per treated unit}, keyed by unit number, with groups separated by {cmd:||}.{p_end}

{marker description}{...}
{title:Description}

{p 4 4 2}
{cmd:multisynth} extends {helpb synth2} from a single treated unit to {bf:multiple} treated units.
It runs the synthetic control method once per treated unit -- each with its own predictor list,
treatment time (which may differ across units, i.e. {bf:staggered adoption}), and optional
matching/evaluation windows -- and then {bf:aggregates} the per-unit treatment effects.

{p 4 4 2}
Effects are aligned on {bf:event time} (periods relative to each unit's own treatment), averaged
with equal weight per unit, and reported over a {bf:common balanced window} equal to the
intersection of the per-unit event spans, so that every event time is averaged over all treated
units (no change in composition along the event-time axis). The overall ATT is the mean of the
averaged gap over the post-treatment event times.

{p 4 4 2}
By default the donor pool for every treated unit is the set of {bf:never-treated} units (all
treated units are excluded); this can be overridden with {opt counit()}. Per-unit estimation
reuses the {helpb synth2} engine, so {helpb synth2} and {helpb synth} options not listed below are
passed through unchanged. The commands {helpb synth} (SSC) and {helpb synth2} are required.

{p 4 4 2}
Per-unit output is {bf:compact} by default (a fit line, the nonzero donor weights, and -- if a
placebo test was requested -- an in-space placebo p-value); specify {opt detail} for the full
{helpb synth2} per-unit output. All results are written to an output folder (see
{help multisynth##outputs:Output files}). See Yan and Chen (2023) for the underlying single-unit
method.

{marker required}{...}
{title:Required settings}

{phang}
{opt trunit(numlist)} a list of one or more treated unit numbers, as given in the panel variable
declared by {helpb xtset} {it:panelvar}.

{phang}
{opt trperiod(unit: # [|| unit: # ...])} the treatment time for {bf:each} treated unit, keyed by
unit number and separated by {cmd:||}; e.g. {cmd:trperiod(3: 1989 || 30: 1985)}. Every unit listed
in {opt trunit()} must have an entry. If all units share the same time, adoption is simultaneous;
otherwise it is staggered.

{phang}
{opt predictors(unit: varlist [|| unit: varlist ...])} the predictor (covariate) list for {bf:each}
treated unit, keyed by unit number and separated by {cmd:||}. Every unit must have its own list;
there is no shared list. Lagged outcomes are written inline in {helpb synth2} style as
{it:depvar}{cmd:(}{it:year}{cmd:)}, e.g. {cmd:cigsale(1988)}.

{marker options}{...}
{title:Options}

{dlgtab:Per-unit periods (keyed; optional)}

{pstd}
Each option below takes a window {bf:per treated unit}, keyed by unit number and separated by
{cmd:||}; e.g. {cmd:xperiod(3: 1980(1)1988 || 30: 1976(1)1984)}. If an option is specified, {bf:every}
treated unit must have an entry -- a single absolute window cannot be correct for units treated at
different times. Omit an option to let each unit use its {helpb synth2} default.

{phang}
{opt xperiod(unit: numlist || ...)} periods over which the covariates in {opt predictors()} are averaged.

{phang}
{opt preperiod(unit: numlist || ...)} pretreatment periods entering the fit.

{phang}
{opt postperiod(unit: numlist || ...)} posttreatment periods over which effects are evaluated.

{phang}
{opt mspeperiod(unit: numlist || ...)} periods over which the mean squared prediction error is minimized.

{dlgtab:Donor pool and output}

{phang}
{opth counit(numlist)} overrides the donor pool with an explicit list of control unit numbers. By
default the donor pool is all {bf:never-treated} units (every unit in {opt trunit()} is excluded),
and the same pool is used for every treated unit.

{phang}
{opt saveresults(folder)} the folder (created if needed, in the current working directory) where
all output files are written. The default is {bf:multisynth_results}.

{phang}
{opt detail} prints the full per-unit {helpb synth2} output (covariate-balance and year-by-year
prediction tables) in addition to the compact per-unit block and the aggregate summary. The
default is compact per-unit output.

{phang}
{opt nofigure} do not build any graphs. By default graphs are built silently (no pop-up windows)
and saved as {bf:.gph} files in the output folder.

{dlgtab:Passed through to the per-unit engine}

{phang}
Any other {helpb synth2} option (e.g. {opt nested}, {opt allopt}, {opt customV()}, {opt margin()},
{opt maxiter()}, {opt sigf()}, {opt bound()}, {opt placebo()}, {opt loo}) is passed through to each
per-unit estimation. When {cmd:placebo(unit)} is requested, the in-space placebo p-value (the
fraction of units whose post/pre MSPE ratio is at least as large as the treated unit's) is reported
per unit and stored in {bf:summary.dta}. See {helpb synth2} for these options. Note that
{opt frame()}, {opt symbol()}, and {opt savegraph()} are {bf:not} available -- {cmd:multisynth}
manages frames and graph saving itself.

{marker outputs}{...}
{title:Output files}

{pstd}
All files are written under the {opt saveresults()} folder (default {bf:multisynth_results}):

{phang}{bf:summary.dta} {hline 2} one row per treated unit: {cmd:trunit uname trperiod att pre_rmspe r2 n_donors t_pre t_post mspe_pval} ({cmd:mspe_pval} is filled only when a placebo test was run).{p_end}
{phang}{bf:aggregate.dta} {hline 2} the equal-weighted average path over the balanced event window: {cmd:reltime}, treated average, synthetic average, ATT, and number of units averaged.{p_end}
{phang}{bf:v_matrix.dta} {hline 2} combined covariate balance / V-weights, one block per unit: {cmd:trunit uname covariate Weight Treated Synthetic Control}.{p_end}
{phang}{bf:weights.dta} {hline 2} combined nonzero donor weights: {cmd:trunit uname donor Weight}.{p_end}
{phang}{bf:unit}{it:#}{bf:_data.dta} {hline 2} the full per-unit series (observed, synthetic, gap) for each treated unit.{p_end}
{phang}{bf:unit}{it:#}{bf:_*.gph}, {bf:aggregate_pred.gph}, {bf:aggregate_eff.gph} {hline 2} per-unit and aggregate graphs (unless {opt nofigure}).{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:multisynth} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(n_treated)}}number of treated units{p_end}
{synopt:{cmd:e(n_donors)}}number of donor (control) units{p_end}
{synopt:{cmd:e(staggered)}}1 if treatment times differ across units, 0 otherwise{p_end}
{synopt:{cmd:e(att_overall)}}overall ATT (mean averaged gap over event time >= 0){p_end}
{synopt:{cmd:e(ev_min)}}lower bound of the balanced event window{p_end}
{synopt:{cmd:e(ev_max)}}upper bound of the balanced event window{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:multisynth}{p_end}
{synopt:{cmd:e(depvar)}}name of the outcome variable{p_end}
{synopt:{cmd:e(predictors)}}the {opt predictors()} specification{p_end}
{synopt:{cmd:e(trperiod)}}the {opt trperiod()} specification{p_end}
{synopt:{cmd:e(trunitlist)}}list of treated unit numbers{p_end}
{synopt:{cmd:e(donors)}}list of donor unit numbers{p_end}
{synopt:{cmd:e(resultsfolder)}}output folder{p_end}

{pstd}
Per-unit results from the underlying {helpb synth2} estimation are written to the output files
listed above rather than left in {cmd:e()}.

{marker references}{...}
{title:References}

{phang}
Abadie, A., A. Diamond, and J. Hainmueller. 2010. Synthetic Control Methods for Comparative Case
Studies: Estimating the Effect of California's Tobacco Control Program.
{it:Journal of the American Statistical Association} 105(490): 493-505.

{phang}
Abadie, A., A. Diamond, and J. Hainmueller. 2015. Comparative Politics and the Synthetic Control
Method. {it:American Journal of Political Science} 59(2): 495-510.

{phang}
Cavallo, E., S. Galiani, I. Noy, and J. Pantano. 2013. Catastrophic Natural Disasters and Economic
Growth. {it:Review of Economics and Statistics} 95(5): 1549-1561.

{phang}
Yan, G. and Q. Chen. 2023. synth2: Synthetic Control Method with Placebo Tests, Robustness Test and
Visualization. {it:The Stata Journal} 23(3): 597-624.

{marker author}{...}
{title:Author}

{pstd}
Ivan Zapata{break}
Texas Tech University{break}
izapatad@ttu.edu{break}

{marker alsosee}{...}
{title:Also see}

{phang}Help: {helpb synth2} (SSC), {helpb synth} (SSC).{p_end}
