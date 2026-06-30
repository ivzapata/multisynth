*! multisynth 0.8.3  (based on synth2 2.1.0 by Yan and Chen)  30jun2026
*  Runs the synthetic control method iteratively over multiple treated units, each with its
*  own predictor list, treatment time (staggered allowed), and matching/evaluation windows,
*  saving per-unit results to a folder, then aggregates the per-unit effects on event
*  (relative) time over a common balanced window into an overall ATT, an aggregate dataset,
*  and event-study graphs. Per-unit output is compact by default (full synth2 output under
*  -detail-); combined v_matrix.dta and weights.dta are saved alongside summary.dta.

program multisynth, eclass
	version 16

	syntax anything, TRUnit(numlist integer sort) TRPeriod(string) PREDictors(string) ///
		[ SAVEResults(string) COUnit(string) FRAme(string) noFIGure ///
		  XPeriod(string) PREPeriod(string) POSTPeriod(string) MSPEPeriod(string) ///
		  DETail * ]

	* require an xtset panel (need panelvar to build the donor pool)
	qui xtset
	if "`r(panelvar)'" == "" | "`r(timevar)'" == "" {
		di as err "panel or time variable missing; use {bf:xtset} {it:panelvar} {it:timevar}"
		exit 459
	}
	local panelvar "`r(panelvar)'"
	local timevar  "`r(timevar)'"

	* outcome variable: only depvar is allowed before the comma
	gettoken depvar extra : anything
	if strtrim("`extra'") != "" {
		di as err "specify only the outcome variable {it:depvar} before the comma; put predictors in {bf:predictors()}"
		exit 198
	}

	* parse predictors(): one keyed group per treated unit, separated by "||"
	*   e.g. predictors( 3: lnincome cigsale(1988) || 30: beer age15to24 )
	local praw `"`predictors'"'
	local seen ""
	while `"`praw'"' != "" {
		if strpos(`"`praw'"', "||") {
			local grp = substr(`"`praw'"', 1, strpos(`"`praw'"', "||") - 1)
			local praw = substr(`"`praw'"', strpos(`"`praw'"', "||") + 2, .)
		}
		else {
			local grp `"`praw'"'
			local praw ""
		}
		local grp = strtrim(`"`grp'"')
		if `"`grp'"' == "" continue
		local cpos = strpos(`"`grp'"', ":")
		if `cpos' == 0 {
			di as err `"predictors(): group "`grp'" is missing a "unit:" key"'
			exit 198
		}
		local key   = strtrim(substr(`"`grp'"', 1, `cpos' - 1))
		local plist = strtrim(substr(`"`grp'"', `cpos' + 1, .))
		local ok : list key in trunit
		if !`ok' {
			di as err "predictors(): unit `key' is not in trunit()"
			exit 198
		}
		local dup : list key in seen
		if `dup' {
			di as err "predictors(): unit `key' is specified more than once"
			exit 198
		}
		if "`plist'" == "" {
			di as err "predictors(): empty predictor list for unit `key'"
			exit 198
		}
		local pred_`key' "`plist'"
		local seen "`seen' `key'"
	}
	* every treated unit must have its own list (no fallback)
	local missing : list trunit - seen
	if "`missing'" != "" {
		di as err "predictors(): no predictors specified for unit(s) `missing'"
		exit 198
	}

	* parse trperiod(): polymorphic, mirroring counit()/predictors():
	*   - one unkeyed integer (e.g. trperiod(1989))  -> common time, all units (simultaneous)
	*   - keyed (e.g. trperiod(3: 1989 || 30: 1985)) -> per-unit times; every unit required
	local traw `"`trperiod'"'
	local tseen ""
	local ttimes ""
	local keyed_trp = 0
	if strpos(`"`traw'"', ":") | strpos(`"`traw'"', "||") local keyed_trp = 1

	if !`keyed_trp' {
		* unkeyed: one common treatment time for every treated unit
		local tval = strtrim(`"`traw'"')
		cap confirm integer number `tval'
		if (_rc) | (`: word count `tval'' != 1) {
			di as err "trperiod(): specify one common time (e.g. trperiod(1989)) or key each unit (e.g. trperiod(3: 1989 || 30: 1985))"
			exit 198
		}
		foreach u of numlist `trunit' {
			local trp_`u' "`tval'"
			local ttimes "`ttimes' `tval'"
		}
	}
	else {
		while `"`traw'"' != "" {
			if strpos(`"`traw'"', "||") {
				local grp = substr(`"`traw'"', 1, strpos(`"`traw'"', "||") - 1)
				local traw = substr(`"`traw'"', strpos(`"`traw'"', "||") + 2, .)
			}
			else {
				local grp `"`traw'"'
				local traw ""
			}
			local grp = strtrim(`"`grp'"')
			if "`grp'" == "" continue
			local cpos = strpos("`grp'", ":")
			if `cpos' == 0 {
				di as err `"trperiod(): group "`grp'" needs a "unit:" key, e.g. trperiod(3: 1989 || 30: 1985)"'
				exit 198
			}
			local key  = strtrim(substr("`grp'", 1, `cpos' - 1))
			local tval = strtrim(substr("`grp'", `cpos' + 1, .))
			local ok : list key in trunit
			if !`ok' {
				di as err "trperiod(): unit `key' is not in trunit()"
				exit 198
			}
			local dup : list key in tseen
			if `dup' {
				di as err "trperiod(): unit `key' is specified more than once"
				exit 198
			}
			cap confirm integer number `tval'
			if _rc {
				di as err "trperiod(): treatment time for unit `key' must be an integer"
				exit 198
			}
			local trp_`key' "`tval'"
			local tseen  "`tseen' `key'"
			local ttimes "`ttimes' `tval'"
		}
		local missing : list trunit - tseen
		if "`missing'" != "" {
			di as err "trperiod(): no treatment time specified for unit(s) `missing'"
			exit 198
		}
	}

	* staggered = treated units do not all share the same treatment time
	local utimes : list uniq ttimes
	local nut : word count `utimes'
	local staggered = (`nut' > 1)

	* parse the absolute-period options, each keyed per unit like predictors():
	*   xperiod( 3: 1980(1)1988 || 30: 1976(1)1984 ), preperiod(...), etc.
	* optional, but if an option is supplied every treated unit must have its own
	* entry (no shared window -- a single absolute window cannot be correct for
	* units treated at different times).
	foreach opt in xperiod preperiod postperiod mspeperiod {
		if `"``opt''"' == "" continue
		local oraw `"``opt''"'
		local oseen ""
		while `"`oraw'"' != "" {
			if strpos(`"`oraw'"', "||") {
				local grp = substr(`"`oraw'"', 1, strpos(`"`oraw'"', "||") - 1)
				local oraw = substr(`"`oraw'"', strpos(`"`oraw'"', "||") + 2, .)
			}
			else {
				local grp `"`oraw'"'
				local oraw ""
			}
			local grp = strtrim(`"`grp'"')
			if `"`grp'"' == "" continue
			local cpos = strpos(`"`grp'"', ":")
			if `cpos' == 0 {
				di as err `"`opt'(): group "`grp'" needs a "unit:" key, e.g. `opt'(3: ... || 30: ...)"'
				exit 198
			}
			local key = strtrim(substr(`"`grp'"', 1, `cpos' - 1))
			local val = strtrim(substr(`"`grp'"', `cpos' + 1, .))
			local ok : list key in trunit
			if !`ok' {
				di as err "`opt'(): unit `key' is not in trunit()"
				exit 198
			}
			local dup : list key in oseen
			if `dup' {
				di as err "`opt'(): unit `key' is specified more than once"
				exit 198
			}
			if "`val'" == "" {
				di as err "`opt'(): empty period list for unit `key'"
				exit 198
			}
			local `opt'_`key' "`val'"
			local oseen "`oseen' `key'"
		}
		local missing : list trunit - oseen
		if "`missing'" != "" {
			di as err "`opt'(): no period specified for unit(s) `missing'"
			exit 198
		}
	}

	if "`frame'" != "" di as txt "note: frame() is ignored; multisynth manages one frame per treated unit."

	* results folder, created in Stata's working directory
	if `"`saveresults'"' == "" local saveresults "multisynth_results"
	capture mkdir `"`saveresults'"'
	local unitdir `"`saveresults'/units"'
	local aggdir  `"`saveresults'/aggregate"'
	capture mkdir `"`unitdir'"'
	capture mkdir `"`aggdir'"'
	di as txt _newline "multisynth: results will be saved under " as res `"`saveresults'/"' as txt " (in `c(pwd)')"

	* donor pool(s). counit() is polymorphic, mirroring predictors():
	*   - empty             -> never-treated units (all treated excluded), shared
	*   - unkeyed numlist   -> one shared donor list applied to every treated unit
	*   - keyed (unit: ...) -> a per-unit donor list (every treated unit required)
	* No treated unit may appear in any donor list.
	qui levelsof `panelvar', local(allunits)
	local nevertr : list allunits - trunit

	local keyed_counit = 0
	if strpos(`"`counit'"', ":") | strpos(`"`counit'"', "||") local keyed_counit = 1

	if "`counit'" == "" {
		foreach u of numlist `trunit' {
			local donors_`u' "`nevertr'"
		}
		local donors "`nevertr'"
	}
	else if !`keyed_counit' {
		capture numlist "`counit'", integer sort
		if _rc {
			di as err "counit(): invalid donor numlist: `counit'"
			exit 198
		}
		local counit "`r(numlist)'"
		local bad : list counit - allunits
		if "`bad'" != "" {
			di as err "counit(): unit(s) `bad' not found in `panelvar'"
			exit 198
		}
		local clash : list trunit & counit
		if "`clash'" != "" {
			di as err "counit(): treated unit(s) `clash' cannot be donors"
			exit 198
		}
		foreach u of numlist `trunit' {
			local donors_`u' "`counit'"
		}
		local donors "`counit'"
	}
	else {
		local craw `"`counit'"'
		local cseen ""
		while `"`craw'"' != "" {
			if strpos(`"`craw'"', "||") {
				local grp = substr(`"`craw'"', 1, strpos(`"`craw'"', "||") - 1)
				local craw = substr(`"`craw'"', strpos(`"`craw'"', "||") + 2, .)
			}
			else {
				local grp `"`craw'"'
				local craw ""
			}
			local grp = strtrim(`"`grp'"')
			if `"`grp'"' == "" continue
			local cpos = strpos(`"`grp'"', ":")
			if `cpos' == 0 {
				di as err `"counit(): group "`grp'" needs a "unit:" key, e.g. counit(3: 1 2 4 || 30: 5 6 7)"'
				exit 198
			}
			local key  = strtrim(substr(`"`grp'"', 1, `cpos' - 1))
			local clst = strtrim(substr(`"`grp'"', `cpos' + 1, .))
			local ok : list key in trunit
			if !`ok' {
				di as err "counit(): unit `key' is not in trunit()"
				exit 198
			}
			local dup : list key in cseen
			if `dup' {
				di as err "counit(): unit `key' is specified more than once"
				exit 198
			}
			if "`clst'" == "" {
				di as err "counit(): empty donor list for unit `key'"
				exit 198
			}
			capture numlist "`clst'", integer sort
			if _rc {
				di as err "counit(): invalid donor numlist for unit `key': `clst'"
				exit 198
			}
			local clst "`r(numlist)'"
			local bad : list clst - allunits
			if "`bad'" != "" {
				di as err "counit(): unit(s) `bad' (for treated unit `key') not found in `panelvar'"
				exit 198
			}
			local clash : list trunit & clst
			if "`clash'" != "" {
				di as err "counit(): treated unit(s) `clash' cannot be donors (for treated unit `key')"
				exit 198
			}
			local donors_`key' "`clst'"
			local cseen "`cseen' `key'"
		}
		local missing : list trunit - cseen
		if "`missing'" != "" {
			di as err "counit(): no donor list specified for unit(s) `missing'"
			exit 198
		}
	}

	* every treated unit must end up with a non-empty donor pool
	foreach u of numlist `trunit' {
		if "`donors_`u''" == "" {
			di as err "no donor units available for treated unit `u'"
			exit 198
		}
	}

	* shared/default expose a single donor count; keyed exposes none
	if `keyed_counit' {
		local donors `"`counit'"'
		local ndonors = .
	}
	else {
		local ndonors : word count `donors'
	}

	*-----------------------------------------------------------------------------
	* run banner
	*-----------------------------------------------------------------------------
	local nt : word count `trunit'
	local timingtxt = cond(`staggered', "staggered", "simultaneous")
	di as txt _newline "{hline 78}"
	di as res " multisynth" as txt "  {c |}  synthetic control for multiple treated units"
	di as txt "{hline 78}"
	di as txt " Outcome          : " as res "`depvar'"
	if `keyed_counit' {
		di as txt " Treated units    : " as res "`nt'" as txt "   Donor pool: " as res "per treated unit (keyed)"
	}
	else {
		local _ptype = cond("`counit'" == "", "never-treated", "shared, user-specified")
		di as txt " Treated units    : " as res "`nt'" as txt "   Donor pool: " as res "`ndonors'" as txt " `_ptype' unit(s)"
	}
	di as txt " Treatment timing : " as res "`timingtxt'"
	di as txt " Output folder    : " as res "`saveresults'/"
	if "`detail'" != "" di as txt " Detail mode      : " as res "on" as txt "  (full per-unit synth2 output)"

	* graph handling: build graphs silently and save as .gph (no pop-up windows)
	if "`figure'" == "" {
		local g0 = c(graphics)
		set graphics off
	}

	* collectors: per-unit summary, combined covariate balance, combined donor weights
	tempname pf
	postfile `pf' int(trunit) str32(uname) int(trperiod) double(att pre_rmspe r2) ///
		int(n_donors t_pre t_post) double(mspe_pval) ///
		using `"`aggdir'/summary.dta"', replace
	tempfile balstack wgtstack
	local balbuilt 0
	local wgtbuilt 0

	* iterate over treated units
	local i 0
	foreach u of numlist `trunit' {
		local ++i
		local uname : label (`panelvar') `u'
		local nd_u : word count `donors_`u''

		* this unit's absolute-period options (keyed per unit; empty if not given)
		local xpopt    ""
		local prepopt  ""
		local postpopt ""
		local mspeopt  ""
		if `"`xperiod'"'    != "" local xpopt    "xperiod(`xperiod_`u'')"
		if `"`preperiod'"'  != "" local prepopt  "preperiod(`preperiod_`u'')"
		if `"`postperiod'"' != "" local postpopt "postperiod(`postperiod_`u'')"
		if `"`mspeperiod'"' != "" local mspeopt  "mspeperiod(`mspeperiod_`u'')"

		* per-unit section header
		di as txt _newline "{hline 78}"
		di as txt " [`i'/`nt'] " as res "`uname'" as txt "   (treated " as res "`trp_`u''" as txt ")"

		* run the engine: quietly by default, fully verbose under -detail-
		capture frame drop _ms_`u'
		if "`detail'" != "" {
			di as txt "{hline 78}"
			_multisynth_unit `depvar' `pred_`u'', trunit(`u') trperiod(`trp_`u'') ///
				counit(`donors_`u'') frame(_ms_`u') `figure' ///
				`xpopt' `prepopt' `postpopt' `mspeopt' `options'
		}
		else {
			quietly _multisynth_unit `depvar' `pred_`u'', trunit(`u') trperiod(`trp_`u'') ///
				counit(`donors_`u'') frame(_ms_`u') `figure' ///
				`xpopt' `prepopt' `postpopt' `mspeopt' `options'
		}

		* save this unit's graphs (.gph) -- built silently via set graphics off
		if "`figure'" == "" {
			foreach g in `e(graph)' {
				capture graph save `g' `"`unitdir'/unit`u'_`g'.gph"', replace
			}
		}
		* save this unit's data (.dta)
		capture frame _ms_`u': save `"`unitdir'/unit`u'_data.dta"', replace

		* in-space placebo p-value (MSPE ratio), only if placebo(unit) was requested
		local pval = .
		capture confirm matrix e(mspe)
		if _rc == 0 {
			if rowsof(e(mspe)) > 1 ///
				mata: st_local("pval", strofreal(mean(st_matrix("e(mspe)")[1,3] :<= st_matrix("e(mspe)")[.,3])))
		}

		* summary row
		post `pf' (`u') ("`uname'") (`trp_`u'') (e(att)) (e(rmse)) (e(r2)) ///
			(`nd_u') (e(T0)) (e(T1)) (`pval')

		* accumulate combined covariate balance (rows = covariates)
		capture confirm matrix e(bal)
		if _rc == 0 {
			tempname B bf
			matrix `B' = e(bal)
			frame create `bf'
			frame `bf' {
				qui set obs `=rowsof(`B')'
				qui svmat double `B', names(col)
				qui gen str32 covariate = ""
				local rn : rownames `B'
				local j 0
				foreach c of local rn {
					local ++j
					qui replace covariate = "`c'" in `j'
				}
				qui gen int trunit = `u'
				qui gen str32 uname = "`uname'"
				order trunit uname covariate
				if `balbuilt' append using `"`balstack'"'
				qui save `"`balstack'"', replace
			}
			frame drop `bf'
			local balbuilt 1
		}

		* accumulate combined donor weights (rows = nonzero donors)
		capture confirm matrix e(U_wt)
		if _rc == 0 {
			tempname W wf
			matrix `W' = e(U_wt)
			frame create `wf'
			frame `wf' {
				qui set obs `=rowsof(`W')'
				qui svmat double `W', names(col)
				qui gen str32 donor = ""
				local rn : rownames `W'
				local j 0
				foreach c of local rn {
					local ++j
					qui replace donor = "`c'" in `j'
				}
				qui gen int trunit = `u'
				qui gen str32 uname = "`uname'"
				order trunit uname donor
				if `wgtbuilt' append using `"`wgtstack'"'
				qui save `"`wgtstack'"', replace
			}
			frame drop `wf'
			local wgtbuilt 1
		}

		* compact per-unit block (under -detail- the engine already printed everything)
		if "`detail'" == "" {
			di as txt "   Control units " as res "`nd_u'" ///
				as txt "    Pre-RMSPE " as res %6.4f e(rmse) ///
				as txt "    R-sq " as res %5.3f e(r2) ///
				as txt "    ATT " as res %8.3f e(att)
			if `pval' < . di as txt "   In-space placebo p-value (MSPE ratio): " as res %5.3f `pval'
			di as txt "   Donor weights:"
			tempname Uc
			matrix `Uc' = e(U_wt)
			local wnames : rownames `Uc'
			mata: synth2_print(tokens("`wnames'")', ("Unit", "U.weight"), st_matrix("`Uc'"), ., 0, 10, 0, 0, 0)
		}

		capture frame drop _ms_`u'
	}
	postclose `pf'

	* restore graphics setting
	if "`figure'" == "" set graphics `g0'

	* save combined balance / weights tables
	if `balbuilt' {
		tempname bsf
		frame create `bsf'
		frame `bsf' {
			use `"`balstack'"', clear
			label variable trunit "treated unit id"
			qui save `"`aggdir'/v_matrix.dta"', replace
		}
		frame drop `bsf'
	}
	if `wgtbuilt' {
		tempname wsf
		frame create `wsf'
		frame `wsf' {
			use `"`wgtstack'"', clear
			label variable trunit "treated unit id"
			qui save `"`aggdir'/weights.dta"', replace
		}
		frame drop `wsf'
	}

	*-----------------------------------------------------------------------------
	* per-unit summary table (synth2 house style; p-value column only if placebo)
	*-----------------------------------------------------------------------------
	tempname sf
	frame create `sf'
	frame `sf' {
		use `"`aggdir'/summary.dta"', clear
		qui count if !missing(mspe_pval)
		local haspval = (r(N) > 0)
		di as txt _newline "{hline 78}"
		di as txt " Per-unit summary" as res "   (saved to `aggdir'/summary.dta)"
		if `haspval' {
			di as txt "{hline 13}{c TT}{hline 54}"
			di as txt %-13s " Unit" "{c |}" %7s "Treat" %8s "Donors" %12s "Pre-RMSPE" %7s "R-sq" %11s "ATT" %9s "p-val"
			di as txt "{hline 13}{c +}{hline 54}"
		}
		else {
			di as txt "{hline 13}{c TT}{hline 45}"
			di as txt %-13s " Unit" "{c |}" %7s "Treat" %8s "Donors" %12s "Pre-RMSPE" %7s "R-sq" %11s "ATT"
			di as txt "{hline 13}{c +}{hline 45}"
		}
		forvalues r = 1/`=_N' {
			local nm = abbrev(uname[`r'], 12)
			if `haspval' {
				di as txt %-13s " `nm'" "{c |}" as res %7.0f trperiod[`r'] %8.0f n_donors[`r'] %12.4f pre_rmspe[`r'] %7.3f r2[`r'] %11.3f att[`r'] %9.3f mspe_pval[`r']
			}
			else {
				di as txt %-13s " `nm'" "{c |}" as res %7.0f trperiod[`r'] %8.0f n_donors[`r'] %12.4f pre_rmspe[`r'] %7.3f r2[`r'] %11.3f att[`r']
			}
		}
		if `haspval' di as txt "{hline 13}{c BT}{hline 54}"
		else         di as txt "{hline 13}{c BT}{hline 45}"
	}
	frame drop `sf'

	*-----------------------------------------------------------------------------
	* Cross-unit aggregation on EVENT (relative) time.
	* Each treated unit's effect path is aligned on reltime = year - trperiod_u,
	* then averaged with equal weight per unit. Averaging is restricted to a
	* COMMON BALANCED window [evlo, evhi] = the intersection of the per-unit event
	* spans, so EVERY event time is averaged over ALL treated units (no change in
	* composition along the event-time axis). The per-unit files keep their full
	* series; only this aggregate is trimmed to the balanced window.
	*-----------------------------------------------------------------------------
	local sym "_"   /* matches the engine's fixed variable-name symbol */
	tempfile stack
	local built 0
	local evlo .
	local evhi .
	foreach u of numlist `trunit' {
		tempname af
		frame create `af'
		frame `af' {
			use `"`unitdir'/unit`u'_data.dta"', clear
			keep if `panelvar' == `u'
			keep `timevar' `depvar' pred`sym'`depvar' tr`sym'`depvar'
			gen double reltime = `timevar' - `trp_`u''
			qui summarize reltime if !missing(tr`sym'`depvar')
			if `built' {
				if `r(min)' > `evlo' local evlo = `r(min)'
				if `r(max)' < `evhi' local evhi = `r(max)'
				append using `"`stack'"'
			}
			else {
				local evlo = `r(min)'
				local evhi = `r(max)'
			}
			save `"`stack'"', replace
		}
		frame drop `af'
		local built 1
	}

	* the balanced window must contain at least one pre- and one post-period
	if `evlo' > `evhi' {
		di as err "cannot aggregate: the treated units' event-time spans do not overlap"
		exit 198
	}
	if `evhi' < 0 {
		di as err "cannot aggregate: the common event window [`evlo', `evhi'] has no post-treatment period"
		exit 198
	}
	if `evlo' > -1 {
		di as err "cannot aggregate: the common event window [`evlo', `evhi'] has no pre-treatment period"
		exit 198
	}

	*-----------------------------------------------------------------------------
	* Pooled permutation inference (Cavallo-Galiani-Noy style).
	* Runs automatically when placebo(unit) was requested -- it needs every
	* donor's gap series, which the placebo pass wrote into unit#_data.dta.
	* For each treated unit we standardize every unit's gap by ITS OWN
	* pre-treatment RMSPE, optionally prune donors whose pre-fit exceeds
	* cutoff x the treated unit's (the same cutoff passed to placebo()), align
	* on event time, then enumerate the full null by crossing one donor per
	* case. Two-sided p = (1 + #{|placebo avg| >= |observed avg|}) / (1 + N).
	*-----------------------------------------------------------------------------
	local doinf  = 0
	local cutoff = .
	local _opts = lower(`"`options'"')
	if strpos(`"`_opts'"', "placebo(") & strpos(`"`_opts'"', "unit") local doinf = 1
	if strpos(`"`_opts'"', "cutoff(") {
		local aft = substr(`"`_opts'"', strpos(`"`_opts'"', "cutoff(") + 7, .)
		local cutoff = real(substr(`"`aft'"', 1, strpos(`"`aft'"', ")") - 1))
	}

	local p_ov2 = .
	local p_ovr = .
	local p_ovl = .
	local n_perms   = .
	if `doinf' {
	quietly {
		local M = `nt'
		* per case: standardize every unit by its own pre-RMSPE, prune by cutoff
		local m 0
		foreach u of numlist `trunit' {
			local ++m
			tempfile case_`m'
			tempname icf
			frame create `icf'
			frame `icf' {
				use `"`unitdir'/unit`u'_data.dta"', clear
				keep `panelvar' `timevar' tr`sym'`depvar'
				drop if missing(tr`sym'`depvar')
				gen double _g2 = (tr`sym'`depvar')^2
				bysort `panelvar': egen double _pm = mean(_g2) if `timevar' < `trp_`u''
				bysort `panelvar': egen double pre_mspe = max(_pm)
				drop _g2 _pm
				gen double std_effect = tr`sym'`depvar' / sqrt(pre_mspe)
				gen int    lead = `timevar' - `trp_`u''
				gen byte   istreated = (`panelvar' == `u')
				if `cutoff' < . {
					qui summarize pre_mspe if istreated, meanonly
					local _trpm = r(mean)
					drop if !istreated & pre_mspe > `cutoff' * `_trpm'
				}
				qui count if !istreated
				if r(N) == 0 {
					noisily di as err "pooled inference: no donors survive cutoff(`cutoff') for unit `u'"
					exit 198
				}
				keep `panelvar' lead istreated std_effect tr`sym'`depvar'
				save `"`case_`m''"', replace
			}
			frame drop `icf'
		}

		* stack the M treated-unit (observed) rows once
		tempfile allreal
		forvalues k = 1/`M' {
			tempname rf
			frame create `rf'
			frame `rf' {
				use `"`case_`k''"', clear
				keep if istreated
				keep lead std_effect tr`sym'`depvar'
				if `k' > 1 append using `"`allreal'"'
				save `"`allreal'"', replace
			}
			frame drop `rf'
		}

		* per-lead pooled test over post event times 0..evhi -> pvalues.dta
		tempname pvf
		postfile `pvf' int(reltime) double(avg_effect p_two p_right p_left n_perms) ///
			using `"`aggdir'/pvalues.dta"', replace
		forvalues L = 0/`evhi' {
			tempname obf
			frame create `obf'
			frame `obf' {
				use `"`allreal'"', clear
				keep if lead == `L'
				collapse (mean) ostd = std_effect (mean) oraw = tr`sym'`depvar'
				local obs_std = ostd[1]
				local obs_abs = abs(ostd[1])
				local obs_raw = oraw[1]
			}
			frame drop `obf'
			forvalues k = 1/`M' {
				tempfile col_`k'
				tempname clf
				frame create `clf'
				frame `clf' {
					use `"`case_`k''"', clear
					keep if !istreated & lead == `L'
					keep std_effect
					rename std_effect e`k'
					save `"`col_`k''"', replace
				}
				frame drop `clf'
			}
			tempname crf
			frame create `crf'
			frame `crf' {
				use `"`col_1'"', clear
				forvalues k = 2/`M' {
					cross using `"`col_`k''"'
				}
				egen double avgp = rowmean(e*)
				local Np = _N
				qui count if abs(avgp) >= `obs_abs'
				local p2 = (r(N) + 1) / (`Np' + 1)
				qui count if avgp >= `obs_std'
				local pr = (r(N) + 1) / (`Np' + 1)
				qui count if avgp <= `obs_std'
				local pl = (r(N) + 1) / (`Np' + 1)
			}
			frame drop `crf'
			post `pvf' (`L') (`obs_raw') (`p2') (`pr') (`pl') (`Np')
		}
		postclose `pvf'
		local n_perms = `Np'   /* identical across leads */

		* windowed (overall ATT) test over reltime >= 0 within the balanced window
		tempname owf
		frame create `owf'
		frame `owf' {
			use `"`allreal'"', clear
			keep if lead >= 0 & lead <= `evhi'
			collapse (mean) cstd = std_effect, by(lead)
			collapse (mean) wstd = cstd
			local obs_w_std = wstd[1]
			local obs_w_abs = abs(wstd[1])
		}
		frame drop `owf'
		forvalues k = 1/`M' {
			tempfile wide_`k'
			tempname wkf
			frame create `wkf'
			frame `wkf' {
				use `"`case_`k''"', clear
				keep if !istreated & lead >= 0 & lead <= `evhi'
				keep `panelvar' lead std_effect
				reshape wide std_effect, i(`panelvar') j(lead)
				forvalues L = 0/`evhi' {
					capture rename std_effect`L' s`k'_`L'
				}
				drop `panelvar'
				save `"`wide_`k''"', replace
			}
			frame drop `wkf'
		}
		tempname wcf
		frame create `wcf'
		frame `wcf' {
			use `"`wide_1'"', clear
			forvalues k = 2/`M' {
				cross using `"`wide_`k''"'
			}
			gen double _wsum = 0
			forvalues L = 0/`evhi' {
				gen double _la = 0
				forvalues k = 1/`M' {
					replace _la = _la + s`k'_`L'
				}
				replace _wsum = _wsum + _la / `M'
				drop _la
			}
			gen double _wavg = _wsum / (`evhi' + 1)
			local Nw = _N
			qui count if abs(_wavg) >= `obs_w_abs'
			local p_ov2 = (r(N) + 1) / (`Nw' + 1)
			qui count if _wavg >= `obs_w_std'
			local p_ovr = (r(N) + 1) / (`Nw' + 1)
			qui count if _wavg <= `obs_w_std'
			local p_ovl = (r(N) + 1) / (`Nw' + 1)
		}
		frame drop `wcf'

		* pooled p-value graphs by event time (styled like synth2's pval*_pboUnit)
		if "`figure'" == "" {
			tempname gpf
			frame create `gpf'
			frame `gpf' {
				use `"`aggdir'/pvalues.dta"', clear
				twoway connected p_two reltime, ///
				       ytitle("two-sided pooled p-values of treatment effects on `depvar'") ///
				       xtitle("event time (relative to treatment)") ///
				       yline(0.05 0.1, lp(dot) lc(black)) ylabel(0(0.1)1) ///
				       title("Pooled Permutation Inference: two-sided p-values") ///
				       name(pooled_pvalTwo_pboUnit, replace) nodraw
				graph save pooled_pvalTwo_pboUnit `"`aggdir'/pooled_pvalTwo_pboUnit.gph"', replace
				twoway connected p_right reltime, ///
				       ytitle("right-sided pooled p-values of treatment effects on `depvar'") ///
				       xtitle("event time (relative to treatment)") ///
				       yline(0.05 0.1, lp(dot) lc(black)) ylabel(0(0.1)1) ///
				       title("Pooled Permutation Inference: right-sided p-values") ///
				       name(pooled_pvalRight_pboUnit, replace) nodraw
				graph save pooled_pvalRight_pboUnit `"`aggdir'/pooled_pvalRight_pboUnit.gph"', replace
				twoway connected p_left reltime, ///
				       ytitle("left-sided pooled p-values of treatment effects on `depvar'") ///
				       xtitle("event time (relative to treatment)") ///
				       yline(0.05 0.1, lp(dot) lc(black)) ylabel(0(0.1)1) ///
				       title("Pooled Permutation Inference: left-sided p-values") ///
				       name(pooled_pvalLeft_pboUnit, replace) nodraw
				graph save pooled_pvalLeft_pboUnit `"`aggdir'/pooled_pvalLeft_pboUnit.gph"', replace
			}
			frame drop `gpf'
		}
	}
	}

	* collapse to the equal-weighted average path over the balanced window
	tempname aggf
	frame create `aggf'
	frame `aggf' {
		use `"`stack'"', clear
		qui keep if reltime >= `evlo' & reltime <= `evhi'
		gen byte inunit = !missing(tr`sym'`depvar')
		collapse (mean) `depvar' pred`sym'`depvar' tr`sym'`depvar' ///
			(sum) n_units = inunit, by(reltime)
		sort reltime
		label variable reltime          "event time (years relative to treatment)"
		label variable `depvar'         "treated average (observed)"
		label variable pred`sym'`depvar' "synthetic average"
		label variable tr`sym'`depvar'  "ATT (average gap)"
		label variable n_units          "number of treated units averaged"
		save `"`aggdir'/aggregate.dta"', replace

		* overall ATT = mean of the average gap over post-treatment event times
		qui summarize tr`sym'`depvar' if reltime >= 0, meanonly
		local att_overall = r(mean)

		* pooled (aggregate) pre-treatment RMSPE: pre-period fit of the averaged synthetic
		qui gen double _pg2 = (tr`sym'`depvar')^2 if reltime < 0
		qui summarize _pg2 if reltime < 0, meanonly
		local pooled_rmspe = sqrt(r(mean))
		qui drop _pg2

		* aggregate event-study graphs (built silently via nodraw, saved as .gph)
		if "`figure'" == "" {
			twoway (line `depvar' reltime, lp(solid)) ///
			       (line pred`sym'`depvar' reltime, lp(dash)), ///
			       title("Average Actual and Synthetic Outcomes") ///
			       xline(-1, lp(dot) lc(black)) ///
			       xtitle("event time (relative to treatment)") ytitle("`depvar'") ///
			       legend(order(1 "Treated avg" 2 "Synthetic avg")) ///
			       name(agg_pred, replace) nodraw
			graph save agg_pred `"`aggdir'/aggregate_pred.gph"', replace
			line tr`sym'`depvar' reltime, ///
			       title("Average Treatment Effect (event time)") ///
			       xline(-1, lp(dot) lc(black)) yline(0, lp(dot) lc(black)) ///
			       xtitle("event time (relative to treatment)") ///
			       ytitle("ATT on `depvar'") legend(off) ///
			       name(agg_eff, replace) nodraw
			graph save agg_eff `"`aggdir'/aggregate_eff.gph"', replace
		}

		* on-screen aggregate event-study table (post-treatment event times)
		qui summarize `depvar' if reltime >= 0, meanonly
		local tmean = r(mean)
		qui summarize pred`sym'`depvar' if reltime >= 0, meanonly
		local smean = r(mean)

		di as txt _newline "{hline 78}"
		di as txt " Average treatment effect on the treated" as res "  (equal-weighted, event time)"
		di as txt "{hline 78}"
		di as txt " Balanced window: " as res "[`evlo', `evhi']" ///
			as txt "    Units averaged: " as res "`nt'" as txt " at every event time"
		di as txt "{hline 11}{c TT}{hline 35}"
		di as txt %11s "Event time" "{c |}" %9s "Treated" %11s "Synthetic" %10s "ATT" %5s "N"
		di as txt "{hline 11}{c +}{hline 35}"
		forvalues r = 1/`=_N' {
			if reltime[`r'] < 0 continue
			di as txt %11.0f reltime[`r'] "{c |}" as res %9.2f `depvar'[`r'] ///
				%11.2f pred`sym'`depvar'[`r'] %10.3f tr`sym'`depvar'[`r'] %5.0f n_units[`r']
		}
		di as txt "{hline 11}{c +}{hline 35}"
		di as txt %11s "Mean" "{c |}" as res %9.2f `tmean' %11.2f `smean' %10.3f `att_overall'
		di as txt "{hline 11}{c BT}{hline 35}"
		di as txt " Pooled pre-treatment RMSPE = " as res %6.4f `pooled_rmspe'
		di as txt "{p 0 6 2}{txt}Note: The average treatment effect over the posttreatment period is{res} " %4.3f `att_overall' "{txt}.{p_end}"
		if !`doinf' di as txt " (pooled p-value: add placebo(unit))"
	}
	frame drop `aggf'

	*-----------------------------------------------------------------------------
	* Pooled inference results (on screen)
	*-----------------------------------------------------------------------------
	if `doinf' {
		local cutmsg = cond(`cutoff' < ., "  Cutoff set at `cutoff'.", "  No cutoff applied.")
		tempname ipf
		frame create `ipf'
		frame `ipf' {
			use `"`aggdir'/pvalues.dta"', clear
			di as txt _newline "{hline 78}"
			di as txt " Pooled permutation inference" as res "   (saved to `aggdir'/pvalues.dta)"
			di as txt "{hline 78}"
			di as txt " P-values compare the observed effect against " as res "`n_perms'" ///
				as txt " placebo combinations." as txt "`cutmsg'"
			di as txt "{hline 12}{c TT}{hline 45}"
			di as txt %12s "Event time" "{c |}" %12s "Avg effect" %11s "p(two)" %11s "p(right)" %11s "p(left)"
			di as txt "{hline 12}{c +}{hline 45}"
			forvalues r = 1/`=_N' {
				di as txt %12.0f reltime[`r'] "{c |}" as res %12.3f avg_effect[`r'] ///
					%11.3f p_two[`r'] %11.3f p_right[`r'] %11.3f p_left[`r']
			}
			di as txt "{hline 12}{c BT}{hline 45}"
			di as txt " p-value of the overall ATT:  p(two) " as res %5.3f `p_ov2' ///
				as txt "   p(right) " as res %5.3f `p_ovr' as txt "   p(left) " as res %5.3f `p_ovl'
			di as txt ""
			di as txt "{p 0 4 2}Note: (1) The two-sided pooled p-value at an event time is the frequency that the absolute placebo average standardized effect is greater than or equal to the absolute observed average standardized effect.{p_end}"
			di as txt "{p 6 4 2}(2) The right-sided (left-sided) pooled p-value is the frequency that the placebo average effect is greater (smaller) than or equal to the observed average effect.{p_end}"
			di as txt "{p 6 4 2}(3) If the estimated average treatment effect is positive, then the right-sided p-value is recommended; whereas the left-sided p-value is recommended if the effect is negative.{p_end}"
		}
		frame drop `ipf'
	}

	* note saved outputs
	di as txt _newline "{hline 78}"
	di as txt " Saved in " as res "`saveresults'/" as txt ":"
	di as txt "   units/     : unit#_data.dta, unit#_*.gph  (one set per treated unit)"
	di as txt "   aggregate/ : summary.dta, aggregate.dta, aggregate_pred.gph, aggregate_eff.gph, v_matrix.dta, weights.dta"
	if `doinf' di as txt "                pvalues.dta" cond("`figure'"=="", ", pooled_pval{Two,Right,Left}_pboUnit.gph", "")

	* minimal stored results
	ereturn clear
	ereturn local cmd           "multisynth"
	ereturn local depvar        "`depvar'"
	ereturn local predictors    `"`predictors'"'
	ereturn local trperiod      `"`trperiod'"'
	ereturn local trunitlist    "`trunit'"
	ereturn local donors        "`donors'"
	ereturn local resultsfolder "`saveresults'"
	ereturn scalar n_treated   = `nt'
	ereturn scalar n_donors    = `ndonors'
	ereturn scalar staggered   = `staggered'
	ereturn scalar att_overall = `att_overall'
	ereturn scalar ev_min      = `evlo'
	ereturn scalar ev_max      = `evhi'
	ereturn scalar pre_rmspe_pooled = `pooled_rmspe'
	if `doinf' {
		ereturn scalar pval_overall       = `p_ov2'
		ereturn scalar pval_overall_right = `p_ovr'
		ereturn scalar pval_overall_left  = `p_ovl'
		ereturn scalar n_perms            = `n_perms'
		ereturn local  pinf_window        "[0, `evhi']"
	}

	di as txt _newline "Finished. Ran multisynth for `nt' treated unit(s); output in `saveresults'/ (units/, aggregate/)."
end

program _multisynth_unit, eclass sortpreserve
	version 16
	preserve
	qui xtset
    if "`r(panelvar)'" == "" | "`r(timevar)'" == "" {
		di as err "panel variable or time variable missing, please use -{bf:xtset} {it:panelvar} {it:timevar}"
		exit 198
    }
	syntax anything, TRUnit(integer) TRPeriod(integer) ///
		[COUnit(numlist min = 1 int sort) ///
		CTRLUnit(numlist min = 1 int sort) ///
		PREPeriod(numlist min = 1 int sort) ///
		POSTPeriod(numlist min = 1 int sort) ///
		XPeriod(passthru) ///
		MSPEPeriod(passthru) ///
		CUStomV(passthru) ///
		nested ///
		allopt ///
		margin(passthru) ///
		maxiter(passthru) ///
		sigf(passthru) ///
		bound(passthru) ///
		placebo(string) ///
		frame(string) ///
		SAVEGraph(string) ///
		noFIGure ///
		]
	local panelVar "`r(panelvar)'"
	local timeVar "`r(timevar)'"
	cap synth
    if _rc== 199 {
	    di as err `"{bf:synth} must be installed (use Stata command "{bf:ssc install synth, replace}")."'
		exit 198
	}
	if "`ctrlunit'" != ""{
		if "`counit'" == ""{
			di as err `"The option {bf:ctrlunit()} is obsolete and replaced by the option {bf:counit()}, but continues to work just like the current option {bf:counit()}."'
			local counit "`ctrlunit'"
		}
		else{
			di as err `"The options {bf:ctrlunit()} and {bf:counit()} can not be specified together, and the current option {bf:counit()} is recommended."'
			exit 198
		}
	}
	/* Check frame */
	if "`frame'" == "" tempname frame
	else {
		capture frame drop `frame'
		qui pwf
		if "`frame'" == "`r(currentframe)'" {
			di as err "invalid frame() -- current frame can not be specified"
			exit 198
		}
		local framename "`frame'"
	}
	/* Check trunit */
	qui levelsof `panelVar', local(unit_n)
	loc check: list trunit in unit_n
	if `check' == 0 {
		di as err "invalid trunit() -- treatment unit not found in {it:panelvar}"
		exit 198
	}
	/* Check counit */
	if "`counit'" != "" {
		loc check: list counit in unit_n
		if `check' == 0 {
			di as err "invalid counit() -- at least one control unit not found in {it:panelvar}"
			exit 198
		}
		loc check: list trunit in counit
		if `check' == 1 {
			di as err "invalid counit() -- treatment unit appears among control units"
			exit 198
		}
		foreach i in `unit_n'{
			loc check: list i in counit
			if `check' == 0 & `i' != `trunit' qui drop if `panelVar' == `i'
		}
	}
	/* Check trperiod */
	qui mata: synth2_levelsof("`timeVar'", "time_n")
	loc check: list trperiod in time_n
	if `check' == 0 {
		di as err "invalid trperiod() -- treatment period not found in {it:timelvar}"
		exit 198
	}
	/* Check preperiod */
	if "`preperiod'" != "" {
		qui mata: synth2_levelsofsel("`timeVar'", "time_pre", st_data(., "`timeVar'"):<`trperiod')
		loc check: list preperiod in time_pre
		if `check' == 0 {
			di as err "invalid preperiod() -- at least one of pretreatment periods that not found in {it:timevar} or not ahead of treatment period"
			exit 198
		}
		foreach i in `time_pre'{
			loc check: list i in preperiod
			if `check' == 0 qui drop if `timeVar' == `i'
		}
	}
	/* Check postperiod */
	if "`postperiod'" != "" {
		qui mata: synth2_levelsofsel("`timeVar'", "time_post", st_data(., "`timeVar'") :>= `trperiod')
		loc check: list posof "`trperiod'" in postperiod
		if `check' != 1 {
			di as err "invalid postperiod() -- treatment period should be the first period of posttreatment periods"
			exit 198
		}
		loc check: list postperiod in time_post
		if `check' == 0 {
			di as err "invalid postperiod() -- at least one of posttreatment periods not found in {it:timevar}"
			exit 198
		}
		foreach i in `time_post'{
			loc check: list i in postperiod
			if `check' == 0 qui drop if `timeVar' == `i'
		}
	}
	/* Symbol used in generated variable names (fixed to "_") */
	local sym "_"
	gettoken depvar indepvars : anything
	local indepvars = strltrim("`indepvars'")
	qui ds
	mata: synth2_abstract("`anything'", "`r(varlist)'", "`depvar'")
	frame put `panelVar' `timeVar' `depvar' `covariates', into(`frame')
	local graphlist ""
	frame `frame'{
		tempvar panelVarStr timeVarStr
		{
			capture decode `panelVar', gen(`panelVarStr')
			if _rc  qui tostring `panelVar', gen(`panelVarStr') usedisplayformat force
			else qui replace `panelVarStr' = subinstr(`panelVarStr', " ", "", .)
			qui replace `panelVarStr' = strtoname(`panelVarStr', 0) 
		}
		qui levelsof `panelVarStr', local(unit_all) clean
		qui levelsof `panelVarStr' if `panelVar' != `trunit', local(unit_ctrl) clean
		qui levelsof `panelVarStr' if `panelVar' == `trunit', local(unit_tr) clean
		
		qui tostring `timeVar', gen(`timeVarStr') usedisplayformat force
		mata: synth2_slevelsof("`timeVarStr'", "time_all")
		mata: synth2_slevelsofsel("`timeVarStr'", "time_pre", st_data(., "`timeVar'") :< `trperiod')
		mata: synth2_slevelsofsel("`timeVarStr'", "time_post", st_data(., "`timeVar'") :>= `trperiod')
		mata: synth2_slevelsofsel("`timeVarStr'", "time_tr", st_data(., "`timeVar'") :== `trperiod')
	}
	qui cap varabbrev synth `anything', trunit(`trunit') trperiod(`trperiod') `xperiod' `mspeperiod' `customV' `margin' `maxiter' `sigf' `bound' `nested' `allopt'
	if (_rc){
		error _rc
		exit
	}
	matrix weight_vars = vecdiag(e(V_matrix))'
	matrix V_wt = e(V_matrix)
	matrix colname weight_vars = Weight
	frame `frame'{
		capture gen pred`sym'`depvar' = .
		label variable pred`sym'`depvar' "prediction of `depvar'
		mata: synth2_insertMatrix("`panelVar'", "`timeVar'", `trunit', st_matrix("e(Y_synthetic)"), "pred`sym'`depvar'")
		capture gen tr`sym'`depvar' = `depvar' - pred`sym'`depvar'
		label variable tr`sym'`depvar' "treatment effect on `depvar'"
		di as txt "Fitting results in the pretreatment periods:"
		mata: synth2_sum("`panelVar'", "`timeVar'", `trunit', `trperiod', "`unit_tr'", "`time_tr'", cols(tokens("`unit_ctrl'")), ///
			cols(tokens("`indepvars'")), st_data(., "pred`sym'`depvar'"), st_data(., "tr`sym'`depvar'"), 0, .)
		matrix balance = e(X_balance), J(rowsof(e(X_balance)), 1 ,.)
		matrix colnames balance = Treated Synthetic Control
	}
	if "`preperiod'" == "" qui mata: synth2_levelsofsel("`timeVar'", "preperiod", st_data(., "`timeVar'") :< `trperiod')
	synth2_balance `indepvars', panelVar(`panelVar') timeVar(`timeVar') trunit(`trunit') `xperiod' preperiod(`preperiod') frame(`frame')
	matrix coljoinbyname balance = weight_vars balance
	mata: temp = st_matrixrowstripe("balance"); ///
		st_matrix("balance", (st_matrix("balance")[., 1..2], st_matrix("balance")[.,3], ((st_matrix("balance")[.,3] :/ st_matrix("balance")[.,2]) :- 1) :* 100, st_matrix("balance")[.,4], ((st_matrix("balance")[.,4] :/ st_matrix("balance")[.,2]) :- 1) :* 100)); ///
		st_matrixrowstripe("balance", temp);
	mata: st_matrixcolstripe("balance", (J(6, 1, ""), ("Vweight", "Treated", "ValueSyntheticControl", "BiasSyntheticControl", "ValueAverageControl", "BiasAverageControl")'))
	tempname balanceFrame
	cap frame create `balanceFrame'
	frame `balanceFrame'{
		qui svmat balance, name(col)
		mata: st_sstore(., st_addvar("strL", "variable"), st_matrixrowstripe("balance")[., 2])
		cap replace BiasSyntheticControl = BiasSyntheticControl
		cap replace BiasAverageControl = BiasAverageControl
		graph dot (asis) BiasAverageControl BiasSyntheticControl, over(variable) yline(0, lcolor(gs10))  marker(1, ///
			msymbol(O)) marker(2, msymbol(X)) ytitle("Standardized % Bias across Covariates") ///
			legend(lab(2 "Synthetic Control") lab(1 "Average Control")) name(bias, replace) title("Covariate Balance") nodraw
			//legend(pos(5) ring(0) lab(1 "Control") lab(2 "Synthetic")) name(bias, replace) nodraw
		local graphlist = "`graphlist' bias"
	}
	di _newline "{p 0 0 2}{txt}Covariate balance in the pretreatment periods:{p_end}"
	mata: synth2_print4(st_matrixrowstripe("balance")[., 2], ("Covariate", "V.weight", "Treated ", "Value     ", "Bias", "Value    ", "Bias "), ///
		(" Synthetic Control", " Average Control"), st_matrix("balance"), ., 0, (8, 10, 10, 7, 10, 7), 0, 0)
	di `"{p 0 6 2}{txt}Note: "V.weight" is the optimal covariate weight in the diagonal of V matrix.{p_end}"'
	di `"{p 6 6 2}{txt}"Synthetic Control" is the weighted average of donor units with optimal weights.{p_end}"'
	di `"{p 6 6 2}{txt}"Average Control" is the simple average of all control units with equal weights.{p_end}"'
	di
	tempname weightVarsFrame
	cap frame create `weightVarsFrame'
	frame `weightVarsFrame'{
		qui svmat weight_vars, name(col)
		mata: st_sstore(., st_addvar("strL", "variable"), st_matrixrowstripe("weight_vars")[., 2])
		qui sum Weight
		qui local ymax= r(max)*1.1
		graph hbar (asis) Weight, over(variable, sort(Weight) descending) ytitle(weight) name(weight_vars, replace) blabel(total, format(%9.4f)) ysc(r(0 `ymax')) title(Optimal Covariate Weights) nodraw
		local graphlist = "`graphlist' weight_vars"
	}
	di "{p 0 6 2}{txt}Optimal Unit Weights:{p_end}"
	frame `frame'{
		mata: synth2_weight("e(W_weights)")
	}
	tempname weightUnitFrame
	cap frame create `weightUnitFrame'
	frame `weightUnitFrame'{
		qui svmat weight_unit, name(col)
		mata: st_sstore(., st_addvar("strL", "unit"), st_matrixrowstripe("weight_unit")[., 2])
		qui sum Weight
		qui local ymax= r(max)*1.1
		graph hbar (asis) Weight, over(unit, sort(Weight) descending) ytitle(weight) name(weight_unit, replace) title(Optimal Unit Weights) blabel(total, format(%9.4f)) ysc(r(0 `ymax')) nodraw 
		local graphlist = "`graphlist' weight_unit"
	}
	di
	frame `frame'{
		di _newline as txt "Prediction results in the posttreatment periods:"
		mata: synth2_print(st_sdata(., "`timeVarStr'"), ("Time", "Actual Outcome", " Synthetic Outcome", " Treatment Effect"), ///
			st_data(., "`depvar' pred`sym'`depvar' tr`sym'`depvar'"), ///
			selectindex((st_data(., "`panelVar'") :== `trunit') :& (st_data(., "`timeVar'") :>= `trperiod')), 1, (13, 17, 16), 0, 0, 1)
		if("`figure'" == ""){
			qui mata: synth2_levelsof("`timeVar'", "temp")
			loc pos: list posof "`trperiod'" in temp
			loc pos = `pos' - 1
			loc xline: word `pos' of `temp'
			twoway (line `depvar' `timeVar') (line pred`sym'`depvar' `timeVar', lpattern(dash)) if `panelVar' == `trunit', ///
					title("Actual and Synthetic Outcomes") xline(`xline', lp(dot) lc(black)) name(pred, replace) ///
					ytitle(`depvar')  legend(order(1 "Actual" 2 "Synthetic")) nodraw
			line tr`sym'`depvar' `timeVar', xline(`xline', lp(dot) lc(black)) yline(0, lp(dot) lc(black)) ///
					title("Treatment Effects") name(eff, replace) ///
					ytitle("treatment effects on `depvar'") nodraw
			local graphlist = "`graphlist' pred eff"
		}
	}
	/* Implement the placebo test using fake unit and/or time */
	if "`placebo'" != "" {
		ereturn local trperiod "`trperiod'"
		ereturn local trunit "`trunit'"
		if(strpos("`placebo'", "unit") != 0 & strpos("`placebo'", "unit(") == 0) local placebo = subinstr("`placebo'", "unit", "unit(.)", .)
		synth2_placebo `anything', trunit(`trunit') trperiod(`trperiod') ///
			panelVar(`panelVar') timeVar(`timeVar') panelVarStr(`panelVarStr') timeVarStr(`timeVarStr') ///
			unit_all(`unit_all') unit_tr(`unit_tr') unit_ctrl(`unit_ctrl') time_tr(`time_tr') time_all(`time_all') ///
			frame(`frame') sym("`sym'") `xperiod' `mspeperiod' `placebo' `figure' `nested' `allopt' `margin' `maxiter' `sigf' `bound'
		local graphlist = "`graphlist' `e(graphlist)'"
		capture mat pval = e(pval)
	}
	/* Display graphs */
	if "`figure'" == "" {
		if "`savegraph'" == "" foreach graph in `graphlist'{
			capture graph display `graph'
		}
		else{
			di
			ereturn local graphlist "`graphlist'"
			synth2_savegraph `savegraph'
		}
	}
	ereturn clear
	capture ereturn matrix pval = pval
	capture if rowsof(mspe) > 1 ereturn matrix mspe = mspe
	ereturn matrix bal = balance
	ereturn matrix U_wt= weight_unit
	ereturn matrix V_wt = V_wt
	capture ereturn local graph "`graphlist'"
	ereturn local frame "`framename'"
	ereturn local time_post "`time_post'"
	ereturn local time_pre "`time_pre'"
	ereturn local time_tr "`time_tr'"
	ereturn local time_all "`time_all'"
	ereturn local unit_ctrl "`unit_ctrl'"
	ereturn local unit_tr "`unit_tr'"
	ereturn local unit_all "`unit_all'"
	ereturn local indepvars "`indepvars'"
	ereturn local depvar "`depvar'"
	ereturn local varlist "`anything'"
	ereturn local timevar "`timeVar'"
	ereturn local panelvar "`panelVar'"
	ereturn scalar N = _N
	ereturn scalar T = wordcount("`time_all'")
	ereturn scalar T0 = wordcount("`time_pre'")
	ereturn scalar T1 = wordcount("`time_post'")
	ereturn scalar K = wordcount("`indepvars'")
	ereturn scalar J = wordcount("`unit_all'")
	ereturn scalar rmse = rmse
	ereturn scalar r2 = r2
	ereturn scalar att = att

	di _newline as txt "Finished."
end

program synth2_savegraph
	version 16
	preserve
	syntax [anything], [asis replace]
	foreach graph in `e(graphlist)'{
		capture graph display `graph'
		graph save `anything'_`graph', `asis' `replace' 
	}
end

program synth2_placebo, eclass sortpreserve
	version 16
	preserve
	local trperiod = `e(trperiod)'
	local trunit = `e(trunit)'
	syntax anything, trunit(numlist) trperiod(numlist) [panelVar(string) timeVar(string) panelVarStr(string) timeVarStr(string) ///
		unit_all(string) unit_tr(string) unit_ctrl(string) ///
		time_all(string) time_tr(string) frame(string) ///
		xperiod(passthru) ///
		mspeperiod(passthru) ///
		customV(passthru) ///
		margin(passthru) ///
		maxiter(passthru) ///
		sigf(passthru) ///
		bound(passthru) ///
		sym(string) ///
		unit(numlist missingokay) period(numlist min = 1 int sort <`trperiod') Cutoff(numlist min = 1 max = 1 >=1) ///
		show(numlist min = 1 max = 1 >=1) noFIGure nested allopt]
	
	gettoken depvar indepvars : anything
	local graphlist ""
	if("`unit'" != "."){
		local unit_pboList ""
		frame `frame'{
			foreach i in `unit'{
				qui levelsof `panelVarStr' if `panelVar' == `i', local(temp) clean
 				local unit_pboList "`unit_pboList' `temp'"
			}			
		}
	}
	else local unit_pboList "`unit_ctrl'"
	
	qui mata: synth2_levelsof("`timeVar'", "temp")
	loc pos: list posof "`trperiod'" in temp
	loc pos = `pos' - 1
	loc xline: word `pos' of `temp'
	if("`unit'" != ""){
		di _newline as txt "Implementing placebo test using fake treatment unit " _continue
		local tsline ""
		local unit_pboSel ""
		local unit_pboRmv ""
		foreach unit_pbo in `unit_pboList'{
			di as res "`unit_pbo'"  as txt "..." _continue
			qui frame `frame': levelsof `panelVar' if `panelVarStr' == "`unit_pbo'", local(pbounit) clean
			qui cap varabbrev synth `anything', trunit(`pbounit') trperiod(`trperiod') `xperiod' `mspeperiod' `customV' `margin' `maxiter' `sigf' `bound' `nested' `allopt'
			if _rc {
			    error _rc
				exit
			}
			frame `frame'{
				capture gen pred`sym'`depvar' = .
				mata: synth2_insertMatrix("`panelVar'", "`timeVar'", `pbounit', st_matrix("e(Y_synthetic)"), "pred`sym'`depvar'")
				qui replace tr`sym'`depvar' = `depvar'- pred`sym'`depvar' if `panelVar' == `pbounit'
				mata: synth2_sum("`panelVar'", "`timeVar'", `pbounit', `trperiod', "`unit_tr'", "`time_tr'", ., ., st_data(., "pred`sym'`depvar'"), ///
					st_data(., "tr`sym'`depvar'"), 1, ("`cutoff'" == "" ? . : strtoreal("`cutoff'")))
				if(`isRmv' == 0){
					local unit_pboSel "`unit_pboSel' `unit_pbo'"
					local line " `line' (tsline tr`sym'`depvar' if `panelVar' == `pbounit', lp(solid) lc(gs8%20)) "
				}
				else local unit_pboRmv "`unit_pboRmv' `unit_pbo'"
			}
		}
		matrix colnames mspe = PreMSPE PostMSPE RatioPostPre RatioTrCtrl
		matrix rownames mspe = `unit_tr' `unit_pboList'
		di _newline _newline as txt "In-space placebo test results using fake treatment units:"
		mata: synth2_print(tokens("`unit_tr' `unit_pboList'")', ("Unit", "Pre MSPE ", " Post MSPE ", "Post/Pre MSPE", "  Pre MSPE of Fake Unit/Pre MSPE of Treated Unit"), st_matrix("mspe"), ., 0, (9, 9, 13, 24), 0, 0, 0)
		if "`unit_pboRmv'" != "" {
			mata: printf(stritrim(sprintf( ///
			"{p 0 6 2}{txt}Note: (1) Using all control units, the probability of obtaining a post/pretreatment MSPE ratio as large as {res}`unit_tr'{txt}'s is{res}%10.4f{txt}.{p_end}\n", mean((st_matrix("mspe")[1, 3] :<= st_matrix("mspe")[. , 3])))))
			mata: printf(stritrim(sprintf( ///
			"{p 6 6 2}{txt}(2) Excluding control units with pretreatment MSPE {res}`cutoff'{txt} times larger than the treated unit, the probability of obtaining a post/pretreatment MSPE ratio as large as {res}`unit_tr'{txt}'s is {res}%10.4f{txt}.{p_end}\n", mean((st_matrix("mspe_cut")[1, 3] :<= st_matrix("mspe_cut")[. , 3])))))
			di "{p 6 6 2}{txt}(3) The pointwise p-values below are computed by excluding control units with pretreatment MSPE {res}`cutoff'{txt} times larger than the treated unit.{p_end}"
			di "{p 6 6 2}{txt}(4) There are total{res}", wordcount("`unit_pboRmv'"), "{txt}units with pretreatment MSPE {res}`cutoff'{txt} times larger than the treated unit, including {res}`unit_pboRmv'{txt}.{p_end}"
		}
		else mata: printf(stritrim(sprintf( ///
			"{p 0 6 2}{txt}Note: The probability of obtaining a post/pretreatment MSPE ratio as large as {res}`unit_tr'{txt}'s is{res}%10.4f{txt}.{p_end}\n", mean((st_matrix("mspe")[1, 3] :<= st_matrix("mspe")[. , 3])))))
		mata: printf("\n{txt}In-space placebo test results using fake treatment units (continued" + ///
			("`cutoff'" == "" ? "" : (", cutoff = {res}`cutoff'")) + "{txt}):\n")
		frame `frame'{
			mata: synth2_placebo(st_data(., "`timeVar'"), st_sdata(., "`panelVarStr'"), st_sdata(., "`timeVarStr'"), "`unit_tr'", "`unit_pboSel'", `trperiod', st_data(., "tr`sym'`depvar'"))
			di "{p 0 6 2}{txt}Note: (1) The two-sided p-value of the treatment effect for a particular period is defined as the frequency that the absolute values of the placebo effects are greater than or equal to the absolute value of treatment effect.{p_end}"
			di "{p 6 6 2}{txt}(2) The right-sided (left-sided) p-value of the treatment effect for a particular period is defined as the frequency that the placebo effects are greater (smaller) than or equal to the treatment effect.{p_end}"
			di "{p 6 6 2}{txt}(3) If the estimated treatment effect is positive, then the right-sided p-value is recommended; whereas the left-sided p-value is recommended if the estimated treatment effect is negative.{p_end}"
			label variable pvalTwo "two-sided p-value of treatment effect generated by 'placebo unit'"
			label variable pvalRight "right-sided p-value of treatment effect generated by 'placebo unit'"
			label variable pvalLeft "left-sided p-value of treatment effect generated by 'placebo unit'"
		}
		if "`figure'" == ""{
			frame `frame'{
				qui local num = wordcount("`unit_pboSel'") + 1
				twoway `line' (tsline tr`sym'`depvar' if `panelVar' == `trunit', lp(solid)) , ///
					xline(`xline', lp(dot) lc(black)) yline(0, lp(dot) lc(black)) ///
					title("In-space Placebo Test") name(eff_pboUnit, replace) ///
					ytitle("treatment/placebo effects on `depvar'") ///
					legend(order(`num' "Treatment Effect" 2 "Placebo Effect")) nodraw
				local graphlist = "`graphlist' eff_pboUnit"
			}
			tempname placeboUnitframe
			
			frame create `placeboUnitframe'
			frame `placeboUnitframe'{
				qui svmat mspe, name(col)
				mata: st_sstore(., st_addvar("strL", "unit"), tokens("`unit_tr' `unit_pboList'")')
				if ("`show'" != "") {
					qui gsort -RatioPostPre
					qui drop if _n >`show'
				}
				graph hbar (asis) RatioPostPre, over(unit, sort(RatioPostPre) descending label(labsize(vsmall))) ///
				ytitle("Ratios of posttreatment MSPE to pretreatment MSPE") ///
				title("In-space Placebo Test") name("ratio_pboUnit", replace) nodraw
			}
			local graphlist = "`graphlist' ratio_pboUnit"
			frame `frame'{
				twoway connected pvalTwo `timeVar' if `panelVar' == `trunit' & `timeVar' >= `trperiod', ///
					ytitle("two-sided p-values of treatment effects on `depvar'") ///
					yline(0.05 0.1, lp(dot) lc(black)) ylabel(0(0.1)1) ///
					title("In-space Placebo Test") name(pvalTwo_pboUnit, replace) nodraw
				twoway connected pvalRight `timeVar' if `panelVar' == `trunit' & `timeVar' >= `trperiod', ///
					ytitle("right-sided p-values of treatment effects on `depvar'") ///
					yline(0.05 0.1, lp(dot) lc(black)) ylabel(0(0.1)1) ///
					title("In-space Placebo Test") name(pvalRight_pboUnit, replace) nodraw
				twoway connected pvalLeft `timeVar' if `panelVar' == `trunit' & `timeVar' >= `trperiod', ///
					ytitle("left-sided p-values of treatment effects on `depvar'") ///
					yline(0.05 0.1, lp(dot) lc(black)) ylabel(0(0.1)1) ///
					title("In-space Placebo Test") name(pvalLeft_pboUnit, replace) nodraw
				local graphlist = "`graphlist' pvalTwo_pboUnit pvalRight_pboUnit pvalLeft_pboUnit"
			}
		}
	}
	if("`period'" != ""){
		di _newline as txt "Implementing placebo test using fake treatment time " _continue
		foreach pboperiod in `period'{
			qui mata: synth2_levelsof("`timeVar'", "time_n")
			loc check: list pboperiod in time_n
			if `check' == 0 {
				di _newline as err "placebo() invalid -- invalid fake period `pboperiod'"
				exit 198
			}
			frame `frame': qui levelsof `timeVarStr' if `timeVar' == `pboperiod', local(time_pbo) clean
			loc pos: list posof "`pboperiod'" in temp
			loc pos = `pos' - 1
			loc xlinePbo: word `pos' of `temp'
			di as res "`time_pbo'" as txt "..." _continue
			qui cap varabbrev synth `anything', trunit(`trunit') trperiod(`pboperiod') `xperiod' `mspeperiod' `customV' `margin' `maxiter' `sigf' `bound' `nested' `allopt'
			if _rc {
			    error _rc
				exit
			}
			frame `frame'{
				loc pboTimeVar = strtoname("`depvar'`sym'`time_pbo'")
				capture gen pred`sym'`pboTimeVar' = .
				mata: synth2_insertMatrix("`panelVar'", "`timeVar'", `trunit', st_matrix("e(Y_synthetic)"), "pred`sym'`pboTimeVar'")
				capture gen tr`sym'`pboTimeVar' = `depvar' - pred`sym'`pboTimeVar'
				label variable pred`sym'`pboTimeVar' "prediction of `depvar' generated by 'placebo period `time_pbo''"
				label variable tr`sym'`pboTimeVar' "treatment effect on `depvar' generated by 'placebo period `time_pbo''"
				if "`figure'" == "" {
					twoway (line `depvar' `timeVar' if `panelVar' == `trunit') ///
						(line pred`sym'`pboTimeVar' `timeVar' if `panelVar' == `trunit', lpattern(dash)), ///
						title("In-time Placebo Test (fake treatment time = `time_pbo')") xline(`xline' `xlinePbo', lp(dot) lc(black)) ///
						name(pred_pboTime`pboperiod', replace) ytitle(`depvar')  ///
						legend(order(1 "Actual" 2 "Synthetic")) nodraw
					local graphlist = "`graphlist' pred_pboTime`pboperiod'"
					line tr`sym'`pboTimeVar' `timeVar' if `panelVar' == `trunit', ///
						xline(`xline' `xlinePbo', lp(dot) lc(black)) yline(0, lp(dot) lc(black)) ///
						title("In-time Placebo Test (fake treatment time = `time_pbo')") name(eff_pboTime`pboperiod', replace) ///
						ytitle("placebo effects on `depvar'") nodraw
					local graphlist = "`graphlist' eff_pboTime`pboperiod'"
				}
			}
		}
		di
		foreach pboperiod in `period'{
			frame `frame': qui levelsof `timeVarStr' if `timeVar' == `pboperiod', local(time_pbo) clean
			loc pboTimeVar = strtoname("`depvar'`sym'`pboperiod'")
			di _newline as txt "In-time placebo test results using fake treatment time " as res "`time_pbo'" as txt":"
			frame `frame': mata: synth2_print(st_sdata(., "`timeVarStr'"), ("Time", "Actual Outcome", " Synthetic Outcome", " Treatment Effect"), ///
				st_data(., "`depvar' pred`sym'`pboTimeVar' tr`sym'`pboTimeVar'"), ///
				selectindex((st_data(., "`panelVar'") :== `trunit') :& (st_data(., "`timeVar'") :>= `pboperiod')), 1, (13, 17, 16), 0, 0, 0)
		}
	}
	capture ereturn matrix pval = pval
	ereturn local graphlist "`graphlist'"
end

**# Balance
program synth2_balance
	version 16
	preserve
	syntax anything, panelVar(string) timeVar(string) trunit(integer) preperiod(numlist) [xperiod(numlist)] frame(string)
	if "`xperiod'" == "" mata: st_local("ifperiod", "(" + invtokens("`timeVar' == " :+ tokens("`preperiod'"), " | ") + ")")
	else mata: st_local("ifperiod", "(" + invtokens("`timeVar' == " :+ tokens("`xperiod'"), " | ") + ")")
	foreach var of local anything{
		local period ""
		local rownumb = rownumb(balance,"`var'")
		if(strpos("`var'","(") > 0){
			local period = substr("`var'", strpos("`var'", "(") + 1, strlen("`var'") - strpos("`var'","(") - 1)
			local period = subinstr("`period'", "&", " ", .)
			qui numlist "`period'"
			local period "`r(numlist)'"
			local var = substr("`var'", 1, strpos("`var'","(") - 1)
		}
		if("`period'" == "") local period "`ifperiod'"
		else mata: st_local("period", "(" + invtokens("`timeVar' == " :+ tokens("`period'"), " | ") + ")")
		cap qui sum `var' if `period' & `panelVar' != `trunit'
		matrix balance[`rownumb', 3] = r(mean)
	}
end

version 16
mata:
	real matrix synth2_uniqrows(real matrix m){
		tmp = J(0, 1, .)
		for(i = 1;i<=rows(m); i++){
			if(i == 1) tmp = tmp\m[i,.]
			else{
				if(sum(tmp:==m[i,.])==0) tmp = tmp\m[i,.]
			}
		}
		return(tmp)
	}
	string matrix synth2_suniqrows(string matrix m){
		tmp = J(0, 1, "")
		for(i = 1;i<=rows(m); i++){
			if(i == 1) tmp = tmp\m[i,.]
			else{
				if(sum(tmp:==m[i,.])==0) tmp = tmp\m[i,.]
			}
		}
		return(tmp)
	}
	void synth2_levelsof(string scalar varname, string scalar localname){
		st_local(localname, "")
		tmp = synth2_uniqrows(st_data(., varname)); 
		for(i=1; i<=rows(tmp);i++) st_local(localname, st_local(localname) + (i==1?"":" ")+ strofreal(tmp[i]))
	}
	void synth2_levelsofsel(string scalar varname, string scalar localname, real matrix selvar){
		st_local(localname, "")
		tmp = synth2_uniqrows(st_data(selectindex(selvar), varname)); 
		for(i=1; i<=rows(tmp);i++) st_local(localname, st_local(localname) + (i==1?"":" ")+ strofreal(tmp[i]))
	}
	void synth2_slevelsof(string scalar varname, string scalar localname){
		st_local(localname, "")
		tmp = synth2_suniqrows(st_sdata(., varname));
		for(i=1; i<=rows(tmp); i++) st_local(localname, st_local(localname) + (i==1?"":" ") + tmp[i])
	}
	void synth2_slevelsofsel(string scalar varname, string scalar localname, real matrix selvar){
		st_local(localname, "")
		tmp = synth2_suniqrows(st_sdata(selectindex(selvar), varname)); 
		for(i=1; i<=rows(tmp);i++) st_local(localname, st_local(localname) + (i==1?"":" ")+ tmp[i])
	}
	void synth2_abstract(string scalar anything, string scalar varlist, string scalar depvar){
		anything = tokens(usubinstr(usubinstr(anything, "(", " ", .), ")", " ", .))
		covariates = ""
		for(i = 1; i <= cols(anything); i++){
			if((sum(anything[i] :==  tokens(varlist)) > 0) & (anything[i] != depvar)) covariates = covariates + " " + anything[i]
		}
		st_local("covariates", invtokens(uniqrows(tokens(covariates)')'))
	}
	void synth2_insertMatrix(string scalar panelVar, string scalar timeVar, real scalar unit_tr, real matrix M, string scalar insertVar){
		real matrix data
		st_view(data, ., invtokens((panelVar, timeVar, invtokens(insertVar))))
		data[selectindex(data[., 1] :== unit_tr), 3..cols(data)] = M
	}
	void synth2_sum(string scalar panelVar, string scalar timeVar, real scalar trunit, real scalar trperiod, string scalar unit_tr, string scalar time_tr, real scalar J, real scalar K, real matrix respo, real matrix effect, real scalar isplacebo, real scalar cut){
		indexPre = selectindex((st_data(., timeVar) :< trperiod) :& (st_data(., panelVar) :== trunit))
		indexPost = selectindex((st_data(., timeVar) :>= trperiod) :& (st_data(., panelVar) :== trunit))
		// MSPE_pre = mean(effect[indexPre, .] :^ 2)
		MSPE_pre =  st_matrix("e(RMSPE)")^2
		MSPE_post = mean(effect[indexPost, .] :^ 2)
		// MSE = MSPE_pre
		MSE = st_matrix("e(RMSPE)")^2
		// RMSE = sqrt(MSE)
		RMSE = st_matrix("e(RMSPE)")
		MAE = mean(abs(effect[indexPre, .]))
		R2 = 1 - sum(effect[indexPre, .] :^ 2)/sum((respo[indexPre, .] :- mean(respo[indexPre, .])):^ 2)
		if(isplacebo == 0){
			wide = 3
			printf("{hline " + strofreal(wide + 77) + "}\n")
			printf(" {txt}%-24uds : {res}%10uds {space "+ strofreal(wide) + "} {txt}%-24uds : {res}%10uds\n", "Treated Unit", abbrev(unit_tr, 10), "Treatment Time", abbrev(time_tr, 10))
			printf("{hline " + strofreal(wide + 77) + "}\n")
			printf(" {txt}%-24uds =  {res}%9.0f {space "+ strofreal(wide) + "} {txt}%-24uds =  {res}%9.5f\n", "Number of Control Units", J, "Root Mean Squared Error", RMSE)
			printf(" {txt}%-24uds =  {res}%9.0f {space "+ strofreal(wide) + "} {txt}%-24uds =  {res}%9.5f\n", "Number of Covariates", K, "R-squared", R2)
			printf("{hline " + strofreal(wide + 77) + "}\n")
			st_numscalar("mse", MSE)
			st_numscalar("mae", MAE)
			st_numscalar("rmse", RMSE)
			st_numscalar("r2", R2)
			st_matrix("mspe", (MSPE_pre, MSPE_post, MSPE_post/MSPE_pre, 1))
			st_matrix("mspe_cut", (MSPE_pre, MSPE_post, MSPE_post/MSPE_pre, 1))
		}
		else{
			ratio = MSPE_pre/st_matrix("mspe")[1, 1]
			st_matrix("mspe", st_matrix("mspe")\(MSPE_pre, MSPE_post, MSPE_post/MSPE_pre, ratio))
			if((cut == .) | (ratio <= cut)) {
				st_local("isRmv", "0")
				st_matrix("mspe_cut", st_matrix("mspe_cut")\(MSPE_pre, MSPE_post, MSPE_post/MSPE_pre, ratio))
			}
			else st_local("isRmv", "1")
		}
	}
	void synth2_weight(string scalar name){
		rownames = st_matrixrowstripe(name)[., 2]
		M = st_matrix(name)[., .]
		delRownames = rownames[selectindex(M[., 2] :== 0), .]'
		rownames = rownames[selectindex(M[., 2] :> 0), .]
		M = M[selectindex(M[., 2] :> 0), .]
		M = sort(((1..rows(M))', M), -3)
		rownames = rownames[M[., 1], .]
		st_local("loounitlist", invtokens(strofreal(M[., 2]')))
		st_local("unit_loolist", invtokens(rownames'))
		M = M[., 3]
		st_matrix("weight_unit", M)
		st_matrixcolstripe("weight_unit", ("", "Weight"))
		st_matrixrowstripe("weight_unit", (J(rows(M), 1, ""), rownames))
		synth2_print(rownames, ("Unit", "U.weight"), M, ., 0, 10, 0, 0, 0)
		if(cols(delRownames) > 0) printf("{p 0 6 2}{txt}Note: The unit {res}" + invtokens(delRownames) ///
			+ "{txt} in the donor pool " + (cols(delRownames) > 1? "get" : "gets") + " a weight of {res}0{txt}.\n")
	}
	void synth2_print(string matrix rownames, string matrix colnames, real matrix M, real matrix indexRow, real scalar isMean, real matrix wideM, real scalar extend, real scalar isInt, real scalar att){
		if (rows(indexRow) != 1){
			rownames = rownames[indexRow, .]
			M = M[indexRow, .]
		}
		wide = max(udstrlen((rownames[1..rows(rownames),1]\colnames[1])))
		printf(sprintf("{hline %g}{c TT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		
		printf(sprintf(" {txt}%%~%guds {c |}", wide), colnames[1])
		for(j = 2; j <= cols(colnames); j++){
			if(j == cols(colnames) & (udstrlen(colnames[j]) > wideM[j - 1] + 2)){
				printf(sprintf("%%%guds\n", wideM[j - 1] + 2), substr(colnames[j], 1, wideM[j - 1]))
				printf(" {space %g} {c |}", wide)
				for(k = 2; k <= cols(colnames); k++){
					if(k != cols(colnames)){
						printf("{space %g}", wideM[k - 1] + 2)
					}
					else printf(sprintf("%%%guds", wideM[k - 1] + 2), substr(colnames[k], wideM[k - 1] +1))
				}
			}
			else printf(sprintf("%%%guds", wideM[j - 1] + 2), colnames[j])
		}
		printf(sprintf("\n{hline %g}{c +}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		for(i = 1; i <= rows(M); i++){
			printf(sprintf(" {txt}%%%guds {c |}{res}", wide), rownames[i])
			for(j = 1; j <= cols(M); j++){
			    if(isInt == 0) printf(sprintf(" %%%g.4f ", wideM[j]), M[i,j])
				else printf(sprintf(" %%%g.0g ", wideM[j]), M[i,j])
			}
			printf("\n")
		}
		if(isMean == 1){
			printf(sprintf("{hline %g}{c +}", wide + 2))
			for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
			printf(sprintf("{hline %g}\n", extend))
			meanM = mean(M[1..rows(M), .])
			printf(sprintf(" {txt}%%~%guds {c |}{res}", wide), "Mean")
			for(j = 1; j <= cols(M); j++){
				printf(sprintf(" %%%g.4f ", wideM[j]), meanM[., j])
			}
			printf("\n")
		}else if(isMean == 2){
			printf(sprintf("{hline %g}{c +}", wide + 2))
			for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
			printf(sprintf("{hline %g}\n", extend))
			sumM = sum(M[1..rows(M), .])
			printf(sprintf(" {txt}%%~%guds {c |}{res}", wide), "Sum")
			for(j = 1; j <= cols(M); j++){
				printf(sprintf(" %%%g.4f ", wideM[j]), sumM[., j])
			}
			printf("\n")
		}
		printf(sprintf("{hline %g}{c BT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		if(isMean == 1){
			printf(stritrim(sprintf("{p 0 6 2}{txt}Note: The average treatment effect over the posttreatment period is{res} %10.4f{txt}.\n", 
				meanM[., 3])))
			if(att == 1) st_numscalar("att", meanM[., 3])
		}
	}
	void synth2_print2(string matrix rownames, string matrix colnames, string scalar colname, real matrix M, real matrix indexRow, real scalar isMean, real matrix wideM, real scalar extend, real scalar isInt){
		if (rows(indexRow) != 1){
			rownames = rownames[indexRow, .]
			M = M[indexRow, .]
		}
		wide = max(udstrlen((rownames[1..rows(rownames),1]\colnames[1])))
		printf(sprintf("{hline %g}{c TT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		
		printf(sprintf(" {txt}%%~%guds {c |}", wide), colnames[1])
		printf(sprintf("{txt}%%%guds", wideM[1] + 2), colnames[2])
		printf(sprintf("{txt}%%~%guds\n", sum(wideM[2..cols(wideM)]) + 2*(cols(wideM)-1)), colname)
		
		printf(sprintf("{space %g}{c |}", wide + 2))
		for(j = 2; j <= cols(colnames); j++){
			if(j == cols(colnames) & (udstrlen(colnames[j]) > wideM[j - 1] + 2)){
				printf(sprintf("%%%guds\n", wideM[j - 1] + 2), substr(colnames[j], 1, wideM[j - 1]))
				printf(" {space %g} {c |}", wide)
				for(k = 2; k <= cols(colnames); k++){
					if(k != cols(colnames)){
						printf("{space %g}", wideM[k - 1] + 2)
					}
					else printf(sprintf("%%%guds", wideM[k - 1] + 2), substr(colnames[k], wideM[k - 1] +1))
				}
			}
			else{
				if(j == 2) printf("{space %g}", wideM[j - 1] + 2)
				else printf(sprintf("%%%guds", wideM[j - 1] + 2), colnames[j])
			}
		}
		printf(sprintf("\n{hline %g}{c +}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		for(i = 1; i <= rows(M); i++){
			printf(sprintf(" {txt}%%%guds {c |}{res}", wide), rownames[i])
			for(j = 1; j <= cols(M); j++){
			    if(isInt == 0) printf(sprintf(" %%%g.4f ", wideM[j]), M[i,j])
				else printf(sprintf(" %%%g.0g ", wideM[j]), M[i,j])
			}
			printf("\n")
		}
		if(isMean == 1){
			printf(sprintf("{hline %g}{c +}", wide + 2))
			for(j = 1; j <= cols(colnames)-1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
			printf(sprintf("{hline %g}\n", extend))
			meanM = mean(M[1..rows(M), .])
			printf(sprintf(" {txt}%%~%guds {c |}{res}", wide), "Mean")
			for(j = 1; j <= cols(M); j++){
				printf(sprintf(" %%%g.4f ", wideM[j]), meanM[., j])
			}
			printf("\n")
		}
		printf(sprintf("{hline %g}{c BT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		if(isMean == 1) 
			printf(stritrim(sprintf("{txt}Note: The average treatment effect over the posttreatment period is{res} %10.4f{txt}.\n", 
				meanM[., 3])))
	}
	void synth2_print3(string matrix rownames, string matrix colnames, string matrix colname, real matrix M, real matrix indexRow, real scalar isMean, real matrix wideM, real scalar extend, real scalar isInt){
		if (rows(indexRow) != 1){
			rownames = rownames[indexRow, .]
			M = M[indexRow, .]
		}
		wide = max(udstrlen((rownames[1..rows(rownames),1]\colnames[1])))
		printf(sprintf("{hline %g}{c TT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		
		printf(sprintf(" {txt}%%~%guds {c |}", wide), colnames[1])
		printf(sprintf("{txt}%%~%guds", sum(wideM[1..2]) + 4), colname[1])
		printf(sprintf("{txt}%%~%guds\n", sum(wideM[3..4]) + 4), colname[2])
		
		printf(sprintf("{space %g}{c |}", wide + 2))
		for(j = 2; j <= cols(colnames); j++){
			if(j == cols(colnames) & (udstrlen(colnames[j]) > wideM[j - 1] + 2)){
				printf(sprintf("%%%guds\n", wideM[j - 1] + 2), substr(colnames[j], 1, wideM[j - 1]))
				printf(" {space %g} {c |}", wide)
				for(k = 2; k <= cols(colnames); k++){
					if(k != cols(colnames)){
						printf("{space %g}", wideM[k - 1] + 2)
					}
					else printf(sprintf("%%%guds", wideM[k - 1] + 2), substr(colnames[k], wideM[k - 1] +1))
				}
			}
			else printf(sprintf("%%%guds", wideM[j - 1] + 2), colnames[j])
		}
		printf(sprintf("\n{hline %g}{c +}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		for(i = 1; i <= rows(M); i++){
			printf(sprintf(" {txt}%%%guds {c |}{res}", wide), rownames[i])
			for(j = 1; j <= cols(M); j++){
			    if(isInt == 0) printf(sprintf(" %%%g.4f ", wideM[j]), M[i,j])
				else printf(sprintf(" %%%g.0g ", wideM[j]), M[i,j])
			}
			printf("\n")
		}
		if(isMean == 1){
			printf(sprintf("{hline %g}{c +}", wide + 2))
			for(j = 1; j <= cols(colnames)-1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
			printf(sprintf("{hline %g}\n", extend))
			meanM = mean(M[1..rows(M), .])
			printf(sprintf(" {txt}%%~%guds {c |}{res}", wide), "Mean")
			for(j = 1; j <= cols(M); j++){
				printf(sprintf(" %%%g.4f ", wideM[j]), meanM[., j])
			}
			printf("\n")
		}
		printf(sprintf("{hline %g}{c BT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		if(isMean == 1) 
			printf(stritrim(sprintf("{txt}Note: The average treatment effect over the posttreatment period is{res} %10.4f{txt}.\n", 
				meanM[., 3])))
	}
	void synth2_print4(string matrix rownames, string matrix colnames, string matrix colname, real matrix M, real matrix indexRow, real scalar isMean, real matrix wideM, real scalar extend, real scalar isInt){
		if (rows(indexRow) != 1){
			rownames = rownames[indexRow, .]
			M = M[indexRow, .]
		}
		wide = max(udstrlen((rownames[1..rows(rownames),1]\colnames[1])))
		printf(sprintf("{hline %g}{c TT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		
		printf(sprintf(" {txt}%%~%guds {c |}", wide), colnames[1])
		printf(sprintf("{txt}%%%guds", wideM[1] + 2), colnames[2])
		printf(sprintf("{txt}%%%guds", wideM[2] + 2), colnames[3])
		printf(sprintf("{txt}%%~%guds", sum(wideM[3..4]) + 4), colname[1])
		printf(sprintf("{txt}%%~%guds\n", sum(wideM[5..6]) + 4), colname[2])
		
		printf(sprintf("{space %g}{c |}", wide + 2))
		for(j = 2; j <= cols(colnames); j++){
			if(j == cols(colnames) & (udstrlen(colnames[j]) > wideM[j - 1] + 2)){
				printf(sprintf("%%%guds\n", wideM[j - 1] + 2), substr(colnames[j], 1, wideM[j - 1]))
				printf(" {space %g} {c |}", wide)
				for(k = 2; k <= cols(colnames); k++){
					if(k != cols(colnames)){
						printf("{space %g}", wideM[k - 1] + 2)
					}
					else printf(sprintf("%%%guds", wideM[k - 1] + 2), substr(colnames[k], wideM[k - 1] +1))
				}
			}
			else{
				if(j == 2 | j == 3) printf("{space %g}", wideM[j - 1] + 2)
				else printf(sprintf("%%%guds", wideM[j - 1] + 2), colnames[j])
			}
		}
		printf(sprintf("\n{hline %g}{c +}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		for(i = 1; i <= rows(M); i++){
			printf(sprintf(" {txt}%%%guds {c |}{res}", wide), rownames[i])
			for(j = 1; j <= cols(M); j++){
			    if(j != 4 & j != 6) printf(sprintf(" %%%g.4f ", wideM[j]), M[i,j])
				else printf(sprintf(" %%%g.2f%%%%", wideM[j]), M[i,j])
			}
			printf("\n")
		}
		if(isMean == 1){
			printf(sprintf("{hline %g}{c +}", wide + 2))
			for(j = 1; j <= cols(colnames)-1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
			printf(sprintf("{hline %g}\n", extend))
			meanM = mean(M[1..rows(M), .])
			printf(sprintf(" {txt}%%~%guds {c |}{res}", wide), "Mean")
			for(j = 1; j <= cols(M); j++){
				printf(sprintf(" %%%g.4f ", wideM[j]), meanM[., j])
			}
			printf("\n")
		}
		printf(sprintf("{hline %g}{c BT}", wide + 2))
		for(j = 1; j <= cols(colnames) - 1; j++) printf(sprintf("{hline %g}", wideM[j] + 2))
		printf(sprintf("{hline %g}\n", extend))
		if(isMean == 1) 
			printf(stritrim(sprintf("{txt}Note: The average treatment effect over the posttreatment period is{res} %10.4f{txt}.\n", 
				meanM[., 3])))
	}
	void synth2_placebo(real matrix timeVar, string matrix panelVarStr, string matrix timeVarStr, string scalar unit_tr, string scalar unit_pboSel, real scalar trperiod, real matrix respo){
		real matrix pvalM
		indexRow = selectindex((panelVarStr :== unit_tr) :& (timeVar :>= trperiod))
		tr_eff = respo[indexRow,]
		unit_list = tokens(unit_pboSel)
		eff = tr_eff
		for(i = 1; i<= cols(unit_list); i++){
			tempIndexRow = selectindex((panelVarStr :== unit_list[i]) :& (timeVar :>= trperiod))
			eff = (eff, respo[tempIndexRow, .])
		}
		pval = J(rows(eff), 3, .)
		for(i = 1; i <= rows(pval); i++) {
			pval[i, 1] = mean((abs(tr_eff[i, 1]) :<= abs(eff[i, .]))')
			pval[i, 2] = mean((tr_eff[i, 1] :<= eff[i, .])')
			pval[i, 3] = mean((tr_eff[i, 1] :>= eff[i, .])')
		}
		synth2_print2(timeVarStr[indexRow, .], ("Time", "Treatment Effect", "Two-sided ", "Right-sided", "Left-sided"), "p-value of Treatment Effect", (tr_eff, pval), ., 0, (16, 11, 11, 11), 0, 0)
		st_matrix("pval", (tr_eff, pval))
		st_matrixcolstripe("pval", (("", "p-value", "p-value", "p-value")',("Tr.Eff.", "two-sided", "right-sided", "left-sided")'))
		st_matrixrowstripe("pval",(J(rows(timeVarStr[indexRow, .]), 1, ""), timeVarStr[indexRow, .]))
		temp = _st_addvar("float", "pvalTwo")
		temp = _st_addvar("float", "pvalRight")
		temp = _st_addvar("float", "pvalLeft")
		st_view(pvalM, ., ("pvalTwo", "pvalRight", "pvalLeft"))
		pvalM[indexRow, .] = pval
	}
end
* 2.1.0 Add the option sign()
* 2.0.1 Fix some bugs
* 2.0.0 Enhance the functionality of saving graphs
* 1.0.0 Address the compatibility issue of varabbreviation
* 0.0.2 Adjust the input of numlist of covariates
* 0.0.1 Fix the issue of parameter transfer in the placebo test
* 0.0.0 Submit the initial version of synth2