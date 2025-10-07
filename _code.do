
*------------------------------------------------------------------------------------------*
* Two-stage model with constructed instruments
* 1. Gravity model regression for agricultural trade – constructs instruments for inputs
* 2. Two-stage regression for manufactured food production
* Created on May 5, 2025 @ Hyungsun Yim
*------------------------------------------------------------------------------------------*

* Install
ssc install ppmlhdfe, replace
ssc install xtivreg2, replace
ssc install ivreg2, replace

* Set directory
pwd


*************************************
*                                   *
* Gravity model: Agricultural trade * 
*                                   *
*************************************

***
* 1) Run regressions for SCTG 01-03
***
eststo clear
use _gravity, clear

forv i = 1/3 {
	
	// drop shared variable names
    capture drop dry_home dry_orig dry_dest ///
		dry_ext_home dry_mild_home dry_ext_orig ///
		dry_mild_orig dry_ext_dest dry_mild_dest ///
		wet_home wet_orig wet_dest ///
		agprod_orig foodprod_dest ///
		gdd_home gdd_orig gdd_dest prec_home prec_orig prec_dest

****** Drought without threshold
    // generate shared variable names
    clonevar dry_home =  ln_dry_sctg`i'_home
    clonevar dry_orig = ln_dry_sctg`i'_orig
    clonevar dry_dest = ln_dry_sctg`i'_dest

    clonevar wet_home  = ln_wet_sctg`i'_home
    clonevar wet_orig = ln_wet_sctg`i'_orig
    clonevar wet_dest = ln_wet_sctg`i'_dest

    clonevar agprod_orig   = ln_gdpi_sctg`i'
    clonevar foodprod_dest = L_ln_gdpj_sctg`i'

    clonevar gdd_home  = ln_gdd_sctg`i'_home
    clonevar gdd_orig = ln_gdd_sctg`i'_orig
    clonevar gdd_dest = ln_gdd_sctg`i'_dest

    clonevar prec_home  = ln_prec_sctg`i'_home
    clonevar prec_orig = ln_prec_sctg`i'_orig
    clonevar prec_dest = ln_prec_sctg`i'_dest

    // estimate ppml
    qui ppmlhdfe sctg`i'tons dry_home dry_orig dry_dest ///
			wet_home wet_orig wet_dest agprod_orig foodprod_dest ///
			mrt_sctg`i'_lntraveltime mrt_sctg`i'_neighbor mrt_sctg`i'_home ///
			gdd_home gdd_orig gdd_dest prec_home prec_orig prec_dest ln_pop, ///
			abs(FEij FEcit FEcjt) cluster(FEij) nolog d
	
	// store results
    eststo model`i'
	
	// generate predicted values
	// later used to construct instruments (I_ijt)
    predict pr_sctg`i', mu
    replace pr_sctg`i' = 0 if pr_sctg`i' == .
	
	// R2* - correlation coefficient between predicted and observed trade flows
    correl sctg`i'tons pr_sctg`i'
    local r2s = r(rho)^2
	
	// store information
    estadd scalar R2s = `r2s'
    estadd scalar R2p = e(r2_p)
    estadd local MRTs  = "Y"
    estadd local FEij  = "Y"
    estadd local FEcit = "Y"
    estadd local FEcjt = "Y"

****** Drought with threshold
    // generate shared variable names
    clonevar dry_ext_home  = ln_dry_sctg`i'_home_ext
    clonevar dry_mild_home = ln_dry_sctg`i'_home_mild
    clonevar dry_ext_orig = ln_dry_sctg`i'_orig_ext
    clonevar dry_mild_orig= ln_dry_sctg`i'_orig_mild
    clonevar dry_ext_dest = ln_dry_sctg`i'_dest_ext
    clonevar dry_mild_dest= ln_dry_sctg`i'_dest_mild

    // estimate ppml
    qui ppmlhdfe sctg`i'tons dry_ext_home dry_mild_home ///
			dry_ext_orig dry_mild_orig dry_ext_dest dry_mild_dest ///
			wet_home wet_orig wet_dest agprod_orig foodprod_dest ///
			mrt_sctg`i'_lntraveltime mrt_sctg`i'_neighbor mrt_sctg`i'_home ///
			gdd_home gdd_orig gdd_dest prec_home prec_orig prec_dest ln_pop, ///
			abs(FEij FEcit FEcjt) cluster(FEij) nolog d

	// store results
    eststo model`i'_th

	// generate predicted values
	// later used to construct instruments (I_ijt)
    predict pr_sctg`i'_th, mu
    replace pr_sctg`i'_th = 0 if pr_sctg`i'_th == .
	
	// R2*
    correl sctg`i'tons pr_sctg`i'_th
    local r2s_th = r(rho)^2
	
	// store information
    estadd scalar R2s = `r2s_th'
    estadd scalar R2p = e(r2_p)
    estadd local MRTs  = "Y"
    estadd local FEij  = "Y"
    estadd local FEcit = "Y"
    estadd local FEcjt = "Y"
}

// export results
esttab model1 model1_th model2 model2_th model3 model3_th ///
    using "_Results_Gravity.csv", replace ///
    b(3) se(3) star(+ 0.12 * 0.10 ** 0.05 *** 0.01) ///
    label drop(mrt_sctg*) nocons ///
    mtitles("Drought(1)" "Extreme/mild(2)" "Drought(1)" "Extreme/mild(2)" "Drought(1)" "Extreme/mild(2)") ///
	order( dry_home dry_orig dry_dest ///
		   dry_ext_home dry_mild_home dry_ext_orig dry_mild_orig dry_ext_dest dry_mild_dest ///
		   wet_home wet_orig wet_dest agprod_orig foodprod_dest ///
		   gdd_home gdd_orig gdd_dest prec_home prec_orig prec_dest ln_pop) ///
    stats(MRTs FEij FEcit FEcjt R2s R2p N, fmt(%s %s %s %s 3 3 0) ///
          labels("MRTs" "State-pair F.E." "Exporter climate zone-year F.E." "Importer climate zone-year F.E." "R-squared*" "Pseudo R-squared" "N"))

// save predicted trade flows
save _iv_pr, replace



***
* 2) Construct instruments and observed inputs (SCTG 01-03)
***
forv i = 1/3 {	

****** Instruments (I_jt)
    foreach s in "" "_th" {
        local pr = "pr_sctg`i'`s'"
        local sum = "sum_sctg`i'`s'"
        local local = "local_sctg`i'`s'"
        local import = "import_sctg`i'`s'"
		
		// sum all inputs by destination
		bysort dest year: egen `sum' = total(`pr')
		// only orig==dest will be kept
        gen `local' = `pr' 
		// imported inputs will be the rest
		gen `import' = `sum' - `local' 
	}
	
****** Observed inputs (LHS of stage 1)
	bysort dest year: egen sum_sctg`i'_end = total(sctg`i'tons)
	gen local_sctg`i'_end = sctg`i'tons
	gen import_sctg`i'_end = sum_sctg`i'_end - local_sctg`i'_end	
 }

// keep instruments and observed ag inputs
// to be merged with stage 2 dataset (rename dest to state)
keep if orig == dest
ren dest state
keep state year local* import*

save _iv_pr, replace



*******************************************
*                                         *
* Two-stage regression: Manufactured food *
*                                         *                                      
*******************************************

***
* 1) Prepare data
***
use _food_stage2, clear

// merge with instruments and observed inputs
merge 1:1 state year using _iv_pr
keep if _merge == 3 
drop _merge

// normalize inputs by labor for constant returns to scale
forv i = 1/3 {	

****** Instruments
	foreach s in "" "_th" {
		local instCD = "instCD_sctg`i'`s'"
        local local = "local_sctg`i'`s'"
        local import = "import_sctg`i'`s'"
		
		gen `instCD' = log(`local'/labor) + log(`import'/labor)
		}

****** Observed inputs
	gen endCD_sctg`i' = log(local_sctg`i'_end/labor) + log(import_sctg`i'_end/labor)
 }
 

// Intrastate flow of SCTG 02 for Rhode Island (2017) is recorded as zero. 
// We interpolate the single value by taking the average of other years. 
// The value we use is -15.
replace endCD_sctg2 = -14.94 if endCD_sctg2 == .



***
* 2) Run Stage 2
***
encode state, gen(state1)
xtset state1 year

// create state and year dummy variables
tab year, gen(t)
tab state, gen(c)

***** Stage 2
eststo clear
foreach s in "" "_th" {

	xtivreg2 ln_y (endCD_sctg1 endCD_sctg2 endCD_sctg3 = instCD_sctg1`s' instCD_sctg2`s' instCD_sctg3`s') ///
		ln_k_FRB t1-t4, cluster(state) fe first
	
	// labor
	lincom 1 - ln_k_FRB - 2*endCD_sctg1 - 2*endCD_sctg2 - 2*endCD_sctg3
	local labor = r(estimate)
	local laborse = r(se)
	
	// store stage 2 results and information
	eststo stage2`s'
	estadd scalar L = `labor'
	estadd scalar Lse = `laborse'
	estadd local FEj = "Y"
	estadd local FEt = "Y"
	
	// R2* - correlation coefficient between predicted and observed values
	qui ivreg ln_y (endCD_sctg1 endCD_sctg2 endCD_sctg3 = instCD_sctg1`s' instCD_sctg2`s' instCD_sctg3`s') ///
		ln_k_FRB c1-c47 t1-t4, robust
	predict pr`s', xb
	correl pr`s' ln_y	
	local r2s = r(rho)^2
	di `r2s'
	
	// store information
	est restore stage2`s'
	estadd scalar R2s = `r2s'	
}

// label
label var endCD_sctg1 "Animals and fish (SCTG 01)"
label var endCD_sctg2 "Cereal grains (SCTG 02)"
label var endCD_sctg3 "Vegetables fruits (SCTG 03)"
label var ln_k_FRB "Capital"

// export stage 2 results
esttab stage2 stage2_th using "_Results_Stage2.csv", replace ///
    b(3) se(3) star(+ 0.12 * 0.10 ** 0.05 *** 0.01) ///
    label drop(t*) nocons ///
    mtitles("Drought(1)" "Ext/Mild Drought(2)") ///
	order(endCD_sctg1 endCD_sctg2 endCD_sctg3 ln_k_FRB) ///
	stats(L Lse FEj FEt R2s N, fmt(3 3 %s %s 3 0) labels("Labor" " " "State F.E." "Year F.E." "R-squared*" "N"))
	
	
***** Stage 1
foreach s in "" "_th" {

	// drop shared names
	capture drop sctg1 sctg2 sctg3
	
	// shared names
	clonevar sctg1 = instCD_sctg1`s'
	clonevar sctg2 = instCD_sctg2`s'
	clonevar sctg3 = instCD_sctg3`s'

    // run regression
    qui xtivreg2 ln_y (endCD_sctg1 endCD_sctg2 endCD_sctg3 = sctg1 sctg2 sctg3) ///
		ln_k_FRB t1-t4, cluster(state) fe first savefirst
	
    // store stage 1 results
	local i = 1
	matrix f = e(first) // weak iv test 
	
    foreach var in endCD_sctg1 endCD_sctg2 endCD_sctg3 {
	   local swF = f[8,`i'] // SW F-stat
	   local rfF = f[4,`i'] // reduced-form F-stat
	   
       est resto _xtivreg2_`var'
       eststo `var'`s'
	   
	   estadd scalar swF = `swF'
	   estadd scalar rfF = `rfF'
	   estadd local FEj = "Y"
       estadd local FEt = "Y"
	   
	   local ++i
    }
}

// label
label var sctg1 "SCTG 01"
label var sctg2 "SCTG 02"
label var sctg3 "SCTG 03"
label var ln_k_FRB "Capital"

// export stage 1 results
esttab endCD_sctg1 endCD_sctg1_th endCD_sctg2 endCD_sctg2_th endCD_sctg3 endCD_sctg3_th using "_Results_Stage1.csv", replace ///
    b(3) se(3) star(+ 0.12 * 0.10 ** 0.05 *** 0.01) ///
    label drop(t*) nocons ///
    mtitles("Drought(1)" "Ext/Mild Drought(2)" "Drought(1)" "Ext/Mild Drought(2)" "Drought(1)" "Ext/Mild Drought(2)") ///
	order( sctg1 sctg2 sctg3 ln_k_FRB) ///
	stats(FEj FEt swF rfF N, fmt(%s %s 2 2 0) labels("State F.E." "Year F.E." "Conditional SW F-statistics" "Reduced-form F-statistics" "N"))


