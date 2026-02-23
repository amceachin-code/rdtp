clear all
pause on

glob projectpath "/dfs2/temp/amceachi/cde/"
glob datapath "${projectpath}generated_data/"

*set matsize 11000

sysdir set PLUS "/dfs2/temp/amceachi/cde/ado"
sysdir set PERSONAL "/dfs2/temp/amceachi/cde/ado"
sysdir set OLDPLACE "/dfs2/temp/amceachi/cde/ado"

set more off, perm

*********************************************
*********************************************	
***** First pass Tipping point analysis *****
*********************************************
*********************************************
cap log close
log using "${projectpath}logs/RD_algorithm_para_school_linear_spline.smcl", replace



foreach t in  2007 2008 2009 2010 2011  {
	use stuid cds year sname cst_mth_alg Lcst_mth_ss analyticS using "${datapath}RD_sample.dta",clear
	
	keep if analyticS==1 
	drop analyticS

	levelsof cds if year==`t', local(sch) 
	
	global model="1.cut##c.force"
	global b="_b[1.cut]"
	global name="Parametric Linear Slope, Interaction"
	global t = "_b[1.cut]/_se[1.cut]"
	
	
	loc counter=1
	foreach cds of local sch {
		use stuid cds year sname cst_mth_alg Lcst_mth_ss sped analyticS using "${datapath}RD_sample.dta", clear
		keep if cds=="`cds'" & year==`t'

		foreach n in 1 {
			foreach v in beta_max cut_max r2_max tstat_max pctR_max pctL_max nR_max nL_max {
				qui g `v'`n'=-99
			}
		}
		cap levelsof Lcst_mth_ss if Lcst_mth_ss>=225 & Lcst_mth_ss<=500 & cds=="`cds'" & year==`t', local(cst) 
		scalar drop _all
		
		foreach scale in beta r2 cut tstat pctR pctL {
			foreach max in 1 {
				local `scale'_`max'=-99
			}
		}
					
		foreach score of local cst {
			di `score'
			di "_________" `beta_1'
			qui g cut=(Lcst_mth_ss>=`score') if Lcst_mth_ss!=. & cds=="`cds'" & year==`t' 
				scalar c=`score'
			
			cap drop force
			qui g force=Lcst_mth_ss-c if cds=="`cds'" & year==`t' 
			
			est clear
			
			cap reg cst_mth_alg ${model}  if cds=="`cds'" & year==`t' & force>=-75 & force<=75, vce(cluster force)	
			scalar Ftest=e(F)
			scalar beta_test=_b[1.cut]
			g bob=1 if Ftest==. | beta_test==0
			tab bob
			if bob==1 {
			}
			else if e(r2) > `r2_1' {
				local beta_1=round(${b}, .001) 
				local r2_1 = round(e(r2),.0001)
				local cut_1=`score' 
				local tstat_1=round(${t}, .001)
				local pctL_1=round(_b[1.cut]+_b[_cons], .001)
				local pctR_1=round(_b[_cons], .001)
				count if force>=0 & force<=75
				local nR_1=r(N)
				count if force<0 & force>=-75
				local nL_1=r(N)
			}
			else {
			}
			drop cut bob
			
			di "`score' " `beta_1'
		}			
		
		replace beta_max1=`beta_1'  if cds=="`cds'" & year==`t'
		replace r2_max1=`r2_1'  if cds=="`cds'" & year==`t'
		replace cut_max1=`cut_1'  if cds=="`cds'" & year==`t'
		replace tstat_max1=`tstat_1'  if cds=="`cds'" & year==`t'
		replace pctR_max1=`pctR_1'  if cds=="`cds'" & year==`t'
		replace pctL_max1=`pctL_1'  if cds=="`cds'" & year==`t'
		cap replace nL_max1=`nL_1'  if cds=="`cds'" & year==`t'
		cap replace nR_max1=`nR_1'  if cds=="`cds'" & year==`t'
		di "______`cds'______`t'___`counter'___"
		

	preserve
		collapse (mean) beta* r2* cut_* tstat* pctR_* pctL_* nR_* nL_* (firstnm) sname, by(cds year)
		order cds sname year beta* cut_* r2* tstat* pctR_* pctL_* nR_* nL_*
		
		if `counter'==1 {
			save "${datapath}RD_algorithm_para_school_linear_spline_`t'", replace
		}
		else {
			append using "${datapath}RD_algorithm_para_school_linear_spline_`t'"
			save "${datapath}RD_algorithm_para_school_linear_spline_`t'", replace
		}
	restore
		
	loc ++counter
	}
}
compress

foreach t in 2007 2008 2009 2010 2011 {
	use "${datapath}RD_sample.dta",clear
		keep if analyticS==1
		keep if year==`t'
		drop analyticS
		merge m:1 cds year using "${datapath}RD_algorithm_para_school_linear_spline_`t'"
		keep if _merge==3
		drop _merge
		tempfile year`t'
		save `year`t''
}
use `year2007'
append using `year2008' `year2009' `year2010' `year2011'

compress
save "${datapath}RD_algorithm_para_school_linear_spline.dta", replace
log close
	
	
