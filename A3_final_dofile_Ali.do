clear all
set more off
cd "C:\Users\AHM\OneDrive\York MA\5820 trade\Assignment 3\final_wd"

//Collect all necessary data

/*Data files we have:-
		bilateral trade --> Trade data in 10-year intervals from 1960 - 2000 (max 5 observations for every imp-exp pair) (17207 pairs)
		
		bilateral_trade_data_2010 --> Trade data for every pair for year 2010 (23,553 pairs)
		
		bilateral_trade_annual --> Trade data for every year from 1962-2000 (max 39 observations for every pair) (21471 pairs)
		
		
*/

/*#################################################################################
				Question 1: Collecting Data
##################################################################################*/

*Step 0: Add 2010 trade data
use bilateral_trade_data_2010.dta, clear
rename ccode_exrtner ccode_exporter
rename tradevalue trade
append using bilateral_trade.dta
save trade_data.dta, replace

*Step 1: Adding GDP & Population
import excel GDP_POP_Size_10YR_1960-2010.xlsx, firstr clear

replace Time = "" in 1600
replace Time = "" in 1601
drop if Time == ""

rename GDPconstant2015USNYGDP GDP
rename GDPpercapitaconstant2015US GDPpc
rename PopulationtotalSPPOPTOTL Pop
rename LandareasqkmAGLNDTOTL Size

local vars "GDP GDPpc Pop Size"

foreach x of local vars{
	replace `x' = "" if `x'==".."
	destring `x', replace
} 
rename Time year
destring year, replace
save GDP_Pop_Size_10yr.dta, replace

//Make separate files for Importing & Exporting Countries
local ei "exporter importer"

foreach y of local ei {
	use GDP_Pop_Size_10yr.dta, clear
	rename CountryCode ccode_`y'
	
	foreach x of local vars {
		rename `x' `x'_`y'
	}
	
	save G_P_S_`y'.dta, replace
}

//Merging
use trade_data.dta, clear
foreach y of local ei {
	merge m:1 year ccode_`y' using G_P_S_`y'.dta
	drop _merge
}
save data_merged_1.dta, replace

*Step 2: Adding Distance/Proximity Factors
use dist_cepii.dta, clear
rename iso_o ccode_importer
rename iso_d ccode_exporter 
keep ccode_exporter ccode_importer dist contig
merge 1:m ccode_exporter ccode_importer using data_merged_1.dta
drop _merge
save data_merged_2.dta, replace


*Step 3: Adding Trade-Openness
//Exports
use data_merged_2.dta, clear
collapse(sum) trade (mean) GDP_exporter, by (year ccode_exporter)
rename trade exports
rename ccode_exporter ccode
sort ccode
save exporters.dta, replace
//Imports
use data_merged_2.dta, clear
collapse(sum) trade (mean) GDP_importer, by (year ccode_importer)
rename trade imports
rename ccode_importer ccode
sort  ccode
save importers.dta, replace   

merge 1:1 ccode year using exporters.dta
drop _merge
gen imp_exp = (imports + exports)*1000
gen tradeopenness= (imp_exp/GDP_exporter)
drop if ccode==""
keep ccode year imp_exp tradeopenness
save Trade_openness.dta, replace

//Make separate files for Importing & Exporting Countries
foreach y of local ei {
	use Trade_openness.dta, clear
	rename ccode ccode_`y'
	rename tradeopenness T_`y'
	rename imp_exp TotTrade_`y'
	save To_`y'.dta, replace
}


use data_merged_2.dta, clear
foreach y of local ei {
	merge m:1 year ccode_`y' using To_`y'.dta
	drop _merge
}
save data_merged_3.dta, replace


*Step 4: Add Landlocked
//Make separate files for Importing & Exporting Countries
foreach y of local ei {
	use geo_cepii.dta, clear 
	keep iso3 landlocked
	duplicates drop iso3, force
	rename iso3 ccode_`y'
	rename landlocked LL_`y'
	save LL_`y'.dta, replace
}


use data_merged_3.dta, clear
foreach y of local ei {
	merge m:1 ccode_`y' using LL_`y'.dta
	drop _merge
}
save data_merged_4.dta, replace


/*#################################################################################
				Question 2: Estimating Eqn (4)
##################################################################################*/

//Part 1: Separate for each year
use data_merged_4.dta, clear
rename ccode_importer ccode
collapse (mean) Y=GDPpc_importer size_i=Size_importer pop_i=Pop_importer T_i=T_importer, by (year ccode)

g ln_Y_i = ln(Y)
g ln_size_i = ln(size_i)
g ln_pop_i = ln(pop_i)

sort year
forvalues year=1970(10)2010{
	eststo m_`year': reg ln_Y_i T_i ln_size_i ln_pop_i if year==`year', cluster(ccode)
}

esttab m* using Q2-1.csv, replace noomitted

//Part 2: For Whole Panel
egen cc_num = group(ccode)
xtset cc_num year, delta(10)

eststo treg1: xtreg ln_Y_i T_i ln_size_i ln_pop_i, cluster(ccode) //No FE

eststo treg2: reghdfe ln_Y_i T_i ln_size_i ln_pop_i, absorb (ccode year) cluster(ccode) //Country FE


esttab treg* using Q2-2.csv, replace noomitted

/*#################################################################################
				Question 3: Estimating Eqn (4)
##################################################################################*/
//Part 1: Separate for each year
use data_merged_4.dta, clear
drop if ccode_importer=="" | ccode_exporter=="" | year==.
*Removing where i=j
drop if ccode_exporter==ccode_importer

egen cpair = group(ccode_exporter ccode_importer)

g pair_id = cond(ccode_importer < ccode_exporter, ccode_importer + "-" + ccode_exporter, ccode_exporter + "-" + ccode_importer)

sort pair_id year
egen pair_trade = total(trade), by(pair_id year)

//i = importer | j = exporter
g ln_tau_GDP = ln((pair_trade*1000)/GDP_importer)
g ln_dist_ij = ln(dist)
g ln_area_i = ln(Size_importer)
g ln_area_j = ln(Size_exporter)
g ln_pop_i = ln(Pop_importer)
g ln_pop_j = ln(Pop_exporter)

sort year
forvalues year=1970(10)2010{
	eststo r_`year': reg ln_tau_GDP c.ln_dist_ij##contig c.ln_area_i##contig c.ln_area_j##contig c.ln_pop_i##contig c.ln_pop_j##contig LL_importer##contig LL_exporter##contig if year==`year', cluster(cpair) 
	predict T_hat_`year' if year ==`year', xb
	replace T_hat_`year' = exp(T_hat_`year') if year ==`year'
}

esttab r* using Q3-1.csv, replace noomitted

//Part 2: For Whole Panel
xtset cpair year, delta(10)

eststo preg1: xtreg ln_tau_GDP c.ln_dist_ij##contig#i.year c.ln_area_i##contig#i.year c.ln_area_j##contig#i.year c.ln_pop_i##contig#i.year c.ln_pop_j##contig#i.year LL_importer##contig#i.year LL_exporter##contig#i.year, cluster(cpair)
predict T_hat_panel, xb
replace T_hat_panel = exp(T_hat_panel)

eststo preg2: reghdfe ln_tau_GDP c.ln_dist_ij##contig#i.year c.ln_area_i##contig#i.year c.ln_area_j##contig#i.year c.ln_pop_i##contig#i.year c.ln_pop_j##contig#i.year LL_importer##contig#i.year LL_exporter##contig#i.year, absorb(cpair year) cluster(cpair)


esttab preg* using Q3-2.csv, replace noomitted

drop _est_r_1970 _est_r_1980 _est_r_1990 _est_r_2000 _est_r_2010 _est_preg1 _est_preg2

save final_data.dta, replace
/*#################################################################################
				Question 4: __________
##################################################################################*/
//Part 1: Separate for each year
rename ccode_importer ccode

collapse (sum) T_hat_1970 T_hat_1980 T_hat_1990 T_hat_2000 T_hat_2010 T_hat_panel (mean) Y=GDPpc_importer size_i=Size_importer pop_i=Pop_importer T_i=T_importer, by(ccode year)

g ln_Y_i = ln(Y)
g ln_size_i = ln(size_i)
g ln_pop_i = ln(pop_i)

forvalues year=1970(10)2010{
	eststo f_`year': reg T_i T_hat_`year' if year==`year' & T_hat_`year'!=0, cluster(ccode)
	eststo f_`year'_c: reg T_i T_hat_`year' ln_size_i ln_pop_i if year==`year' & T_hat_`year'!=0, cluster(ccode)
	esttab f_`year'* using Q4-`year'.csv, replace noomitted
	
	/*Not working in loop:-
	
	twoway (scatter T_i T_hat_`year' if year==`year', ms(circle_hollow)) (lfit T_i T_hat_`year' if year==`year'), ///
title("Figure 2: Relationship b/w Actual vs Constructed Trade Share `year'") ytitle("Actual T") xtitle("Constructed T") legend(off)

	graph export using `year'fig.pdf, replace*/
}

twoway (scatter T_i T_hat_1970 if year==1970, ms(circle_hollow)) (lfit T_i T_hat_1970 if year==1970), ///
title("Fig 2: Relationship b/w Actual vs Constructed Trade Share 1970        ") ytitle("Actual T") xtitle("Constructed T") legend(off)
//graph export using fig_year1.pdf, replace

twoway (scatter T_i T_hat_1980 if year==1980, ms(circle_hollow)) (lfit T_i T_hat_1980 if year==1980), ///
title("Fig 2: Relationship b/w Actual vs Constructed Trade Share 1980       ") ytitle("Actual T") xtitle("Constructed T") legend(off) 
//graph export using fig_year2.pdf, replace

twoway (scatter T_i T_hat_1990 if year==1990, ms(circle_hollow)) (lfit T_i T_hat_1990 if year==1990), ///
title("Fig 2: Relationship b/w Actual vs Constructed Trade Share 1990       ") ytitle("Actual T") xtitle("Constructed T") legend(off)
//graph export using fig_year3.pdf, replace

twoway (scatter T_i T_hat_2000 if year==2000, ms(circle_hollow)) (lfit T_i T_hat_2000 if year==2000), ///
title("Fig 2: Relationship b/w Actual vs Constructed Trade Share 2000       ") ytitle("Actual T") xtitle("Constructed T") legend(off)
//graph export using fig_year4.pdf, replace

twoway (scatter T_i T_hat_2010 if year==2010, ms(circle_hollow)) (lfit T_i T_hat_2010 if year==2010), ///
title("Fig 2: Relationship b/w Actual vs Constructed Trade Share 2010       ") ytitle("Actual T") xtitle("Constructed T") legend(off)
//graph export using fig_year5.pdf, replace


//Part 2: For Whole Panel
eststo f_panel: reg T_i T_hat_panel if T_hat_panel!=0, cluster(ccode)
eststo f_panel_c: reg T_i T_hat_panel ln_size_i ln_pop_i if T_hat_panel!=0, cluster(ccode)

drop _est_f_1970 _est_f_1980 _est_f_1990 _est_f_2000 _est_f_2010 _est_f_panel _est_f_1970_c _est_f_1980_c _est_f_1990_c _est_f_2000_c _est_f_2010_c

esttab f_panel* using Q4-P.csv, replace noomitted

twoway (scatter T_i T_hat_panel, ms(circle_hollow)) (lfit T_i T_hat_panel), ///
title("Fig 2: Relationship b/w Actual vs Constructed Trade Share Panel      ") ytitle("Actual T") xtitle("Constructed T") legend(off)

graph export fig2_panel.pdf, replace
/*#################################################################################
				Question 5: Comparing OLS and IV
##################################################################################*/
//Part 1: Separate for each year
forvalues year=1970(10)2010{
	eststo frOLS_`year': reg ln_Y_i T_i ln_size_i ln_pop_i if year==`year', cluster(ccode)
	eststo frIV_`year': ivreg ln_Y_i ln_size_i ln_pop_i (T_i = T_hat_`year') if year==`year', cluster(ccode)
}

//Part 2: For Whole Panel
eststo frOLS_panel: reg ln_Y_i T_i ln_size_i ln_pop_i, cluster(ccode)
eststo frIV_panel: ivreg ln_Y_i ln_size_i ln_pop_i (T_i = T_hat_panel), cluster(ccode)

esttab fr* using Q5.csv, replace noomitted

drop _est_frOLS_1970 _est_frOLS_1980 _est_frOLS_1990 _est_frOLS_2000 _est_frOLS_2010 _est_frOLS_panel _est_frIV_1970 _est_frIV_1980 _est_frIV_1990 _est_frIV_2000 _est_frIV_2010 _est_frIV_panel

/*#################################################################################
				Remove Unnecessary Files
##################################################################################*/
rm final_data.dta
rm LL_importer.dta
rm LL_exporter.dta
rm To_exporter.dta
rm To_importer.dta
rm Trade_openness.dta
rm importers.dta
rm exporters.dta
rm G_P_S_importer.dta
rm G_P_S_exporter.dta
rm data_merged_4.dta
rm data_merged_3.dta
rm data_merged_2.dta
rm data_merged_1.dta
rm GDP_Pop_Size_10yr.dta
rm trade_data.dta