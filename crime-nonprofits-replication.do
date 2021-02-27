**********************************************************************************************************************
**********************************************************************************************************************
**********************************************************************************************************************
*** Replication code for "Community and the Crime Decline: The Causal Effect of Local Nonprofits on Violent Crime"
*** Authors: Patrick Sharkey, Gerard Torrats-Espinosa, and Delaram Takyar
*** Created: 08/29/2017
**********************************************************************************************************************
**********************************************************************************************************************
**********************************************************************************************************************

**********************************************************************************************************************
*** Table 2 : Models of long term change
**********************************************************************************************************************
set more off

use "crime-nonprofits-long-term.dta", clear
local controls_d popdens_d asian_d black_d hispanic_d  other_d lesshs_d college_d fborn_d male1524_d poverty_d unemployed_d manufacturing_d

* Run OLS 
*******************************************
foreach crime in murd viol prop {
reg log_`crime'_r_d all_cml_r_d  all_cml_r  `controls_d' [aw=totpop], robust
estimates store ols_d_all_`crime'
}

* Run 2SLS 
*******************************************
* First stage
reg all_cml_r_d iv_cml_r_d all_cml_r   `controls_d' [aw=totpop] , robust
estimates store iv_d_all_first
test iv_cml_r_d
estadd scalar ftest `r(F)'

* Second stage
foreach crime in murd viol prop {
ivreg2  log_`crime'_r_d  (all_cml_r_d = iv_cml_r_d)   `controls_d'   all_cml_r [aw=totpop], robust 
estimates store iv_d_all_`crime'
}

* Generate Table 2
*******************************************
esttab ols_d_all_murd ols_d_all_viol ols_d_all_prop    iv_d_all_first     iv_d_all_murd iv_d_all_viol iv_d_all_prop using  "table2", replace rtf ///
	mgroups("OLS" "2SLS", pattern(1 0 0 1 0 0 0)) ///
	mtitles("Murder" "Violent" "Property" "Community nonprofits" "Murder" "Violent" "Property" )  ///
	drop( _cons) ///
	order(all_cml_r_d iv_cml_r_d  all_cml_r)  ///
	label  eqlabels(none)  ///
	b(3) p(3)  ///
	star(* 0.05 ** 0.01 *** 0.001) ///
	nolines ///
	collabels(none) ///
	cells("b(fmt(3)star)" "se(fmt(3)par)") ///
	stats(ftest N r2_a    , fmt(%9.3fc %9.0fc %9.3fc )   labels(`"F-test IV"' `"Observations"'  `"Adj. R2"')) ///
	title("Table 2: Long-term change estimates for community nonprofits") ///
    note("* p<0.05, ** p<0.01, *** p<0.001 (two-tailed tests). Heteroskedasticty-robust standard errors in parentheses. All models include 1990 population weights.") ///

**********************************************************************************************************************
*** Table 3: Fixed effects models
**********************************************************************************************************************
set more off
use "crime-nonprofits-panel.dta", clear

xtset place_id year
tsset place_id year
local controls  popdens asian black hispanic  other lesshs college fborn male1524 poverty unemployed manufacturing

* 2SLS first stage
*******************************************
	xtreg  all_cml_r   iv_cml_r `controls' y_*   [aw=popweights],fe cluster(place_id)
	estimates store ivfe_first_cml
	test iv_cml_r
	estadd scalar ftest `r(F)'
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

foreach dv in murd viol prop {

*OLS
*******************************************
	xtreg log_`dv'_r_1ld all_cml_r `controls' y_*   [aw=popweights],fe cluster(place_id) 
	estimates store olsfe_`dv'_cml
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

* 2SLS second stage
*******************************************
	xtivreg2 log_`dv'_r_1ld (all_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_cml
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

* 2SLS reduced form
*******************************************
	xtreg log_`dv'_r_1ld iv_cml_r `controls' y_*   [aw=popweights],fe cluster(place_id) 
	estimates store ivfe_`dv'_red_cml
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

}

* Genearte Table 3 
*******************************************
esttab olsfe_murd_cml  olsfe_viol_cml olsfe_prop_cml ivfe_first_cml ivfe_murd_cml ivfe_viol_cml ivfe_prop_cml  using "table3", replace     rtf   ///
	mgroups("OLS" "2SLS", pattern(1 0 0 1 0 0 0)) ///
	mtitles("Murder" "Violent" "Property" "Community nonprofits" "Murder" "Violent" "Property" )  ///
	drop(y_* _cons) ///
	order(all_cml_r iv_cml_r )  ///
	label  eqlabels(none)  ///
	b(3) p(3)  ///
	star(* 0.05 ** 0.01 *** 0.001) ///
	nolines ///
	collabels(none) ///
	cells("b(fmt(3)star)" "se(fmt(3)par)") ///
	stats(ftest N r2_a  cityfe yearfe , fmt(%9.3fc %9.0fc %9.3fc )   labels(`"F-test IV"' `"Observations"'  `"Adj. R2"'  `"City fixed effects"' `"Year fixed effects"')) ///
	title("Table 3: OLS and IV fixed effects estimates for all community nonprofits") ///
	note("* p < 0.05, ** p < 0.01, *** p < 0.001 (two-tailed tests). Standard errors clustered by city in parentheses. All models include 1990 population weights.")


**********************************************************************************************************************
*** Table 4: Robustness tests
**********************************************************************************************************************
foreach i in murd viol prop {

**  Baseline 
 xtreg log_`i'_r_1ld  all_cml_r `controls' y_*   [aw=popweights],fe cluster(place_id)

	mat ols_base_`i' = J(1,3,.)
	mat ols_base_`i'[1,1] = _b[all_cml_r]
	mat ols_base_`i'[1,2] = _se[all_cml_r]
	mat ols_base_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r) `controls' y_*   [aw=popweights],fe cluster(place_id)

	mat iv_base_`i' = J(1,3,.)
	mat iv_base_`i'[1,1] = _b[all_cml_r]
	mat iv_base_`i'[1,2] = _se[all_cml_r]
	mat iv_base_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))

*** Without population weights
 xtreg log_`i'_r_1ld  all_cml_r `controls' y_*  ,fe cluster(place_id)

	mat ols_nowgt_`i' = J(1,3,.)
	mat ols_nowgt_`i'[1,1] = _b[all_cml_r]
	mat ols_nowgt_`i'[1,2] = _se[all_cml_r]
	mat ols_nowgt_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r) `controls' y_*  ,fe cluster(place_id)

	mat iv_nowgt_`i' = J(1,3,.)
	mat iv_nowgt_`i'[1,1] = _b[all_cml_r]
	mat iv_nowgt_`i'[1,2] = _se[all_cml_r]
	mat iv_nowgt_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))


** Baseline with year by region FE
 xtreg log_`i'_r_1ld  all_cml_r `controls' year_reg_*  [aw=popweights],fe cluster(place_id) 

	mat ols_reg_`i' = J(1,3,.)
	mat ols_reg_`i'[1,1] = _b[all_cml_r]
	mat ols_reg_`i'[1,2] = _se[all_cml_r]
	mat ols_reg_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r) `controls' year_reg_*  [aw=popweights],fe cluster(place_id) 
	
	mat iv_reg_`i' = J(1,3,.)
	mat iv_reg_`i'[1,1] = _b[all_cml_r]
	mat iv_reg_`i'[1,2] = _se[all_cml_r]
	mat iv_reg_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))


** Lagged dependent variable
 xtreg log_`i'_r_1ld log_`i'_r_1lg all_cml_r `controls' y_*   [aw=popweights],fe cluster(place_id)
	
	mat ols_ldv_`i' = J(1,3,.)
	mat ols_ldv_`i'[1,1] = _b[all_cml_r]
	mat ols_ldv_`i'[1,2] = _se[all_cml_r]
	mat ols_ldv_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  log_`i'_r_1lg (all_cml_r=iv_cml_r) `controls' y_*   [aw=popweights],fe cluster(place_id)
	
	mat iv_ldv_`i' = J(1,3,.)
	mat iv_ldv_`i'[1,1] = _b[all_cml_r]
	mat iv_ldv_`i'[1,2] = _se[all_cml_r]
	mat iv_ldv_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))

**  100 largest cities
 xtreg log_`i'_r_1ld  all_cml_r `controls' y_*   [aw=popweights] if poprank<=100,fe cluster(place_id)
	
	mat ols_100_`i' = J(1,3,.)
	mat ols_100_`i'[1,1] = _b[all_cml_r]
	mat ols_100_`i'[1,2] = _se[all_cml_r]
	mat ols_100_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r) `controls' y_*   [aw=popweights]  if poprank<=100,fe cluster(place_id)
	
	mat iv_100_`i' = J(1,3,.)
	mat iv_100_`i'[1,1] = _b[all_cml_r]
	mat iv_100_`i'[1,2] = _se[all_cml_r]
	mat iv_100_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))

** Excluding 100 largest cities
 xtreg log_`i'_r_1ld  all_cml_r `controls' y_*   [aw=popweights] if poprank>100,fe cluster(place_id)
	
	mat ols_exc100_`i' = J(1,3,.)
	mat ols_exc100_`i'[1,1] = _b[all_cml_r]
	mat ols_exc100_`i'[1,2] = _se[all_cml_r]
	mat ols_exc100_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r) `controls' y_*   [aw=popweights]  if poprank>100,fe cluster(place_id)
	
	mat iv_exc100_`i' = J(1,3,.)
	mat iv_exc100_`i'[1,1] = _b[all_cml_r]
	mat iv_exc100_`i'[1,2] = _se[all_cml_r]
	mat iv_exc100_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))

** Baseline with incarceration 
 xtreg log_`i'_r_1ld  all_cml_r  log_incarceraton  `controls' y_*   [aw=popweights],fe cluster(place_id)

	mat ols_inc_`i' = J(1,3,.)
	mat ols_inc_`i'[1,1] = _b[all_cml_r]
	mat ols_inc_`i'[1,2] = _se[all_cml_r]
	mat ols_inc_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r)  log_incarceraton  `controls' y_*   [aw=popweights],fe cluster(place_id)
	
	mat iv_inc_`i' = J(1,3,.)
	mat iv_inc_`i'[1,1] = _b[all_cml_r]
	mat iv_inc_`i'[1,2] = _se[all_cml_r]
	mat iv_inc_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))

**  Baseline with police officers 
 xtreg log_`i'_r_1ld  all_cml_r swornftime_r  `controls' y_*   [aw=popweights],fe cluster(place_id)

	mat ols_pol_`i' = J(1,3,.)
	mat ols_pol_`i'[1,1] = _b[all_cml_r]
	mat ols_pol_`i'[1,2] = _se[all_cml_r]
	mat ols_pol_`i'[1,3] = 2 * ttail(e(df_r), abs(_b[all_cml_r]/_se[all_cml_r]))

 xtivreg2 log_`i'_r_1ld  (all_cml_r=iv_cml_r) swornftime_r  `controls' y_*   [aw=popweights],fe cluster(place_id)
	
	mat iv_pol_`i' = J(1,3,.)
	mat iv_pol_`i'[1,1] = _b[all_cml_r]
	mat iv_pol_`i'[1,2] = _se[all_cml_r]
	mat iv_pol_`i'[1,3] = 2 * ttail(e(Fdf2), abs(_b[all_cml_r]/_se[all_cml_r]))

}

*** Stack results from different models by crime type
foreach i in murd viol prop {
	mat ols_all_`i' = ols_base_`i'
	mat ols_all_`i' = ols_all_`i'\ols_nowgt_`i'
	mat ols_all_`i' = ols_all_`i'\ols_reg_`i'
	mat ols_all_`i' = ols_all_`i'\ols_ldv_`i'
	mat ols_all_`i' = ols_all_`i'\ols_100_`i'
	mat ols_all_`i' = ols_all_`i'\ols_exc100_`i'
	mat ols_all_`i' = ols_all_`i'\ols_inc_`i'
	mat ols_all_`i' = ols_all_`i'\ols_pol_`i'
}

foreach i in murd viol prop {
	mat iv_all_`i' = iv_base_`i'
	mat iv_all_`i' = iv_all_`i'\iv_reg_`i'
	mat iv_all_`i' = iv_all_`i'\iv_nowgt_`i'
	mat iv_all_`i' = iv_all_`i'\iv_ldv_`i'
	mat iv_all_`i' = iv_all_`i'\iv_100_`i'
	mat iv_all_`i' = iv_all_`i'\iv_exc100_`i'
	mat iv_all_`i' = iv_all_`i'\iv_inc_`i'
	mat iv_all_`i' = iv_all_`i'\iv_pol_`i'
}

*** Append horizonytally
local cols = rowsof(ols_all_murd)
mat empty = J(`cols',1,.)

mat mat_ols = ols_all_murd, ols_all_viol,  ols_all_prop
mat mat_iv =  iv_all_murd,  iv_all_viol,  iv_all_prop
mat mat_ols_iv = mat_ols,  mat_iv


mat rownames mat_ols_iv = /*
*/  "(1) Baseline" /*
*/  "(2) Without population weights" /*
*/  "(3) Region-by-year fixed effects"/*
*/  "(4) Lagged dependent variable" /*
*/  "(5) Only 100 largest cities" /*
*/  "(6) Excluding 100 largest cities" /*
*/  "(7) Adding state incarceration" /*
*/  "(8) Adding police (1992-2008)" 

**** Display OLS and 2SLS coefficients from each model in one table
* Subset coeffs and std error columns (don't need the p-values)
local cols = rowsof(ols_all_murd)
foreach mat in ols_all_murd iv_all_murd ols_all_viol iv_all_viol  ols_all_prop iv_all_prop {
mat `mat' = `mat'[1..`cols',1],`mat'[1..`cols',2]  
}

*** Create matrix with star signifficance for each set of results
di invttail(5000,.05/2)
di invttail(5000,.01/2)
di invttail(5000,.001/2)
foreach mat in ols_all_murd iv_all_murd ols_all_viol iv_all_viol  ols_all_prop iv_all_prop {
	local cols = rowsof(`mat')
	matrix stars_`mat' = J(`cols',2,0)
	forvalues k = 1/`cols' {
/* change to matrix stars_`mat'[`k',2] to get stars next to SEs (* .05,. ** .01, *** .001) */	matrix stars_`mat'[`k',1] =  (abs(`mat'[`k',1]/`mat'[`k',2]) >1.96) + (abs(`mat'[`k',1]/`mat'[`k',2]) > 2.58)  + (abs(`mat'[`k',1]/`mat'[`k',2]) > 3.29)
	}
} 

*** Create 6 tables seperately (3 crimes x 2 models) and then merge them)
foreach mat in ols_all_murd iv_all_murd ols_all_viol iv_all_viol  ols_all_prop iv_all_prop {
mat rownames `mat' = /*
*/  "(1) Baseline" /*
*/  "(2) Without population weights" /*
*/  "(3) Region-by-year fixed effects"/*
*/  "(4) Lagged dependent variable" /*
*/  "(5) Only 100 largest cities" /*
*/  "(6) Excluding 100 largest cities" /*
*/  "(7) Adding state incarceration" /*
*/  "(8) Adding police (1992-2008)" 

frmttable, statmat(`mat') substat(1) sdec(3) annotate(stars_`mat') asymbol(*,**,***) varlabels
}

* Genearte Table 4
*******************************************
frmttable using  table4a, statmat(ols_all_murd) substat(1) sdec(3) annotate(stars_ols_all_murd)  asymbol(*,**,***)  title("Table 4: Robustness to alternate model specifications and sample choices")  ctitle("",  "Murder OLS" )    replace
frmttable using  table4a, statmat(iv_all_murd ) substat(1) sdec(3) annotate(stars_iv_all_murd) asymbol(*,**,***) merge ctitle("",   "Murder 2SLS" )    
frmttable using  table4a, statmat(ols_all_viol) substat(1) sdec(3) annotate(stars_ols_all_viol) asymbol(*,**,***) merge ctitle("",  "Violent OLS" )   
frmttable using  table4a, statmat(iv_all_viol) substat(1) sdec(3) annotate(stars_iv_all_viol) asymbol(*,**,***) merge ctitle("",  "Violent 2SLS" )    
frmttable using  table4a, statmat(ols_all_prop) substat(1) sdec(3) annotate(stars_ols_all_prop) asymbol(*,**,***) merge ctitle("",  "Property OLS" )   
frmttable using  table4a,statmat(iv_all_prop) substat(1) sdec(3) annotate(stars_iv_all_prop) asymbol(*,**,***) merge   ctitle("",  "Property 2SLS" ) note("* p < 0.05, ** p < 0.01, *** p < 0.001 (two-tailed tests). Standard errors clustered by city in parentheses. All models include 1990 population weights")


**********************************************************************************************************************
*** Table A1: Lagged effects of crime on IV
**********************************************************************************************************************
areg iv_new_r log_prop_r_4lg y_*   [aw=popweights],a(place_id) cluster(place_id)  
predict residlag, residuals

foreach dv in murd viol prop {

	xtreg iv_new_r log_`dv'_r_4lg y_*   [aw=popweights] if residlag!=.,fe cluster(place_id) 
	estimates store lag_`dv'_nc
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

	xtreg iv_new_r log_`dv'_r_4lg   `controls' y_*   [aw=popweights] if residlag!=.,fe cluster(place_id) 
	estimates store lag_`dv'_c
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"
}
lab var log_murd_r_4lg "Log murder (3-year lag)"
lab var log_viol_r_4lg "Log violent (3-year lag)"
lab var log_prop_r_4lg "Log property (3-year lag)"


* Genearte Table A1
*******************************************
esttab lag_murd_nc lag_murd_c lag_viol_nc lag_viol_c  lag_prop_nc lag_prop_c  using "tableA1", replace rtf ///
	mgroups("Arts, medical, & environmental nonprofits per 100,000", pattern(1 0 0 0 0 0))  ///
	nomtitles  ///
	drop(y_* _cons) ///
	order(log_murd_r_4lg log_viol_r_4lg log_prop_r_4lg )  ///
	label  eqlabels(none)  ///
	b(3) p(3)  ///
	star(* 0.05 ** 0.01 *** 0.001) ///
	nolines ///
    nonotes ///
	collabels(none) ///
	cells("b(fmt(3)star)" "se(fmt(3)par)") ///
	stats(N r2_a  cityfe yearfe , fmt(%9.0fc %9.3fc )   labels(`" Observations"'  `"Adj. R2"'  `"City fixed effects"' `"Year fixed effects"')) ///
	title("Table A1: Lagged effects of crime rates on the instrument") ///
	note("* p < 0.05, ** p < 0.01, *** p < 0.001 (two-tailed tests). Standard errors clustered by city in parentheses. All models include 1990 population weights.")


**********************************************************************************************************************
*** Table A2: Newly added nonprofits
**********************************************************************************************************************
* 2SLS first stage
	xtreg  all_new_r   iv_new_r `controls' y_*   [aw=popweights],fe cluster(place_id)
	estimates store ivfe_first_new
	test iv_new_r
	estadd scalar ftest `r(F)'
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

foreach dv in murd viol prop {
*OLS
	xtreg log_`dv'_r_1ld all_new_r `controls' y_*   [aw=popweights],fe cluster(place_id) 
	estimates store olsfe_`dv'_new
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"

* 2SLS second stage
	xtivreg2 log_`dv'_r_1ld (all_new_r  = iv_new_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_new
	estadd local cityfe "Yes"
	estadd local yearfe "Yes"
}

* Generate Table A2
******************************************
esttab olsfe_murd_new  olsfe_viol_new olsfe_prop_new ivfe_first_new ivfe_murd_new ivfe_viol_new ivfe_prop_new  using "tableA2", replace rtf ///
	mgroups("OLS" "2SLS", pattern(1 0 0 1 0 0 0))  ///
	mtitles("Murder" "Violent" "Property" "Community nonprofits" "Murder" "Violent" "Property" )  ///
	drop(y_* _cons) ///
	order(all_new_r iv_new_r )  ///
	label  eqlabels(none)  ///
	b(3) p(3)  ///
	star(* 0.05 ** 0.01 *** 0.001) ///
	nolines ///
    nonotes ///
	collabels(none) ///
	cells("b(fmt(3)star)" "se(fmt(3)par)") ///
	stats(ftest N r2_a  cityfe yearfe , fmt(%9.2fc %9.0fc %9.3fc %9.0fc)   labels(`"F-test IV"' `" Observations"'  `"Adj. R2"'  `"City fixed effects"' `"Year fixed effects"')) ///
	title("Table A2: OLS and IV fixed effects estimates using newly added nonprofits") ///
	note("* p < 0.05, ** p < 0.01, *** p < 0.001 (two-tailed tests). Standard errors clustered by city in parentheses. All models include 1990 population weights.")


**********************************************************************************************************************
*** Figure 4: Lags and leads
**********************************************************************************************************************
* lag and lead the IV and the treatment
sort place_id year
foreach var in iv_new iv_cml  all_new all_cml {
	* Leads
	by place_id: gen `var'_r_5ld =   `var'_r[_n+5]
	by place_id: gen `var'_r_4ld =   `var'_r[_n+4]
	by place_id: gen `var'_r_3ld =   `var'_r[_n+3]
	by place_id: gen `var'_r_2ld =   `var'_r[_n+2]
	by place_id: gen `var'_r_1ld =   `var'_r[_n+1]
	* Lags
	by place_id: gen `var'_r_5lg =   `var'_r[_n-5]
	by place_id: gen `var'_r_4lg =   `var'_r[_n-4]
	by place_id: gen `var'_r_3lg =   `var'_r[_n-3]
	by place_id: gen `var'_r_2lg =   `var'_r[_n-2]
	by place_id: gen `var'_r_1lg =   `var'_r[_n-1]
}

* 2SLS
foreach l in r_5ld r_4ld r_3ld r_2ld r_1ld r  r_1lg r_2lg r_3lg r_4lg r_5lg {
ren all_new_`l' treat
ren  iv_new_`l'  iv
xtivreg2  log_murd_r   ( treat = iv)  y_*   [aw=popweights],fe cluster(place_id)
estimates store ivmurd_`l'

xtivreg2  log_viol_r      ( treat = iv)   y_*   [aw=popweights],fe cluster(place_id)
estimates store ivviol_`l'

xtivreg2  log_prop_r      ( treat = iv)  y_*   [aw=popweights],fe cluster(place_id)
estimates store ivprop_`l'

ren treat all_new_`l' 
ren iv  iv_new_`l'  
}

* Generate figure 4
******************************************
coefplot    ivviol_r_3lg ivviol_r_2lg ivviol_r_1lg ivviol_r  ivviol_r_1ld ivviol_r_2ld ivviol_r_3ld  ,  ///
	scheme(s1mono)  plotregion(lcolor(none)) ///
	keep(treat) ///
	/// order( )  ///
	vertical ///
	bycoefs ///
	title("", size(medium)) ///
	ytitle("Effect of community nonprofits formation" "on violent crime rate (2SLS)", size(medium)) ///
	ylabel(-.5(.25).25,labsize(medsmall)) ///
	yscale(titlegap(*5)) ///.
	bylabels( " -3y          -2y          -1y           0y          +1y          +2y          +3y ") ///
	xtitle("Lags and leads between nonprofits and crime", size(medium))  ///
	xlabel(,labsize(medsmall)) ///
	xscale(titlegap(*5))  ///
	ciopts(lpatt(solid)lcol(blue)) ///
	color(blue) msymbol(O) msize(1.2) mfcolor(blue) ///
	yline(0,lpattern(dash))   ///
	legend(off) 
	graph export figure4.pdf, replace


**********************************************************************************************************************
*** Figure 5: Effects for each nonprofit type
**********************************************************************************************************************
lab var all_cml_r "All community organizations"
lab var crime_cml_r "Crime prevention"
lab var nhood_cml_r "Neighborhood development"
lab var substance_cml_r "Substance abuse programs"
lab var jobs_cml_r "Workforce development"
lab var youth_cml_r "Youth organizations"

foreach dv in murd viol prop {

* 2SLS second stage
	xtivreg2 log_`dv'_r_1ld (all_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_all

	xtivreg2 log_`dv'_r_1ld (crime_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_crime

	xtivreg2 log_`dv'_r_1ld (nhood_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_nhood

	xtivreg2 log_`dv'_r_1ld (substance_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_substance

	xtivreg2 log_`dv'_r_1ld (jobs_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_jobs

	xtivreg2 log_`dv'_r_1ld (youth_cml_r  = iv_cml_r)   `controls' y_*   [aw=popweights],fe  cluster(place_id) 
	estimates store ivfe_`dv'_youth
}


* Generate figure 5
******************************************
coefplot  ivfe_murd_crime ivfe_murd_nhood   ivfe_murd_substance ivfe_murd_jobs ivfe_murd_youth, bylabel(Murder) ///
	||  ivfe_viol_crime ivfe_viol_nhood   ivfe_viol_substance ivfe_viol_jobs ivfe_viol_youth, bylabel(Violent crime)  ///		
	||  ivfe_prop_crime ivfe_prop_nhood   ivfe_prop_substance ivfe_prop_jobs ivfe_prop_youth, bylabel(Property crime)  ///
	||, drop( `controls' y_* _cons) ///
	scheme(s1mono)  plotregion(lcolor(none)) ///
	subtitle(, size(small)  margin(small)   lstyle(none) fcolor(none)) ///
	nokey ///
	nooffsets ///
	mlabel format(%9.3f) mlabposition(11) mlabsize(2.2) mlabgap(*1.5)  ///		
	xlabel(-.3(.1).05,labsize(small) format(%9.1f)) /// 
	coeflabels(,wrap(15) nobreak) ///
	ylabel(,labsize(small)) ////
	grid(none) ///
	xline(0,lpattern(dash) lcolor(black) lwidth(.15)) ///
	color(blue) msymbol(O) msize(.8) mfcolor(blue)  ///
	ciopts(lpatt(solid)lcol(blue)) ///
	byopts( cols(3)) 
graph export figure5.pdf, replace

