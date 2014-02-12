clear all
set more off

/*This program runs an analysis of CRE property prices in particular U.S. markets.
The main objective is to create summary statistics and measures of dispersion of
distressed assets and intensity of repeat buyers.
*/



********
**The loop begins
********
local cbsas 29820 33100 47900 35620 41860

local d29820 "Las Vegas, NV"
local d33100 "Miami, FL"
local d47900 "Washington, D.C."
local d35620 "New York, NY"
local d41860 "San Francisco, CA"


foreach poi in `cbsas' {
	use if CBSA_cd == `poi' using ~/analysts/tanya/FR1013.dta, clear

	****Keep only SALES

	keep if TransType_tx == "Sale"


	****Make sure we don't have any duplicates
	duplicates drop Deal_id Status_dt, force


	****Date to average by
	gen year = year(Status_dt)


	****Make a numeric Buyer ID
	egen B_id = group(BuyerName1)
	drop if B_id >= .


	****Variable for length of window in terms of days: 3 years here
	local window = 1095


	****Find the repeat buyers
	gen keeper = 1
	egen kp = total(keeper), by(B_id)


	****Create an indicator of repeat buyers
	gen keep = .
	replace keep = 1 if kp > 1
	replace keep = 0 if kp == 1

	gen indic = .
	levelsof B_id if keep == 1, local(b_ids)


	****Find and flag the repeat buyers within our specified window
	foreach b in `b_ids' {	
		levelsof Status_dt if B_id == `b', local(b_times)
		foreach y in `b_times' {
			local plus = 0
			sum B_id if (Status_dt < `y' & Status_dt >= (`y' - `window') & B_id == `b')
			local plus = `plus' + `r(N)'
			replace indic = 1 if (B_id == `b' & Status_dt == `y' & `plus' != 0)
			}
		}

	replace indic = 0 if indic >= .


	****Calculate percentages of repeat buyers by year
	gen not = 1
	egen numer = total(indic), by(year)
	egen denom = total(not), by(year)
	gen stat =  (numer/denom)*100


	****Now we make our graphs of the repeat buyers as a percentage of all buyers
	preserve

	duplicates drop year stat, force
	sort year

	twoway line stat year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
	title(Repeat Buyers as a Percentage of New Sales) subtitle("`d`poi'': 3-year Window")
	graph export ~/analysts/tanya/`poi'_cre_repeat_3y.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_repeat_3y.dta, replace
	
	restore


	****Now we make our graphs for number of new transactions
	drop kp
	egen kp = total(keeper), by(year)

	preserve

	duplicates drop year kp, force
	sort year

	twoway line kp year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Transactions) ///
	title(New Transactions: Number) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_cre_newtran_num.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_newtran_num.dta, replace
	
	restore


	****Now we make our graphs for dollar amount of new transactions
	replace Price = Price/1000000000

	preserve

	egen tot_amount = total(Price), by(year)
	duplicates drop year tot_amount, force
	sort year

	twoway line tot_amount year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle("$ Billions") ///
	title(New Transactions: Dollar Amount) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_cre_newtran_amt.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_newtran_amt.dta, replace

	restore


	****Now we make our graphs for percentage of dollar amount sold to repeat buyers
	drop numer denom stat
	egen numer = total(indic*Price), by(year)
	egen denom = total(not*Price), by(year)
	gen stat =  (numer/denom)*100

	preserve

	duplicates drop year stat, force
	sort year

	twoway line stat year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
	title(Percentage of Dollar Amount Sold to Repeat Buyers) subtitle("`d`poi'': 3-year Window")
	graph export ~/analysts/tanya/`poi'_cre_repeat_3y_perc.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_repeat_3y_perc.dta, replace

	restore


	****Here we flag repeat sales of the same property
	drop indic
	gen indic = .
	egen g_pkid = group(PropertyKey_id)

	levelsof g_pkid, local(p_ids)

	foreach b in `p_ids' {	
		levelsof Status_dt if g_pkid == `b', local(p_times)
		foreach y in `p_times' {
			local plus = 0
			sum g_pkid if (Status_dt < `y' & g_pkid == `b')
			local plus = `plus' + `r(N)'
			replace indic = 1 if (g_pkid == `b' & Status_dt == `y' & `plus' != 0)
			}
		}

	replace indic = 0 if indic >= .


	****Now we graph percentage of repeat sales by property
	drop not numer denom stat
	gen not = 1
	egen numer = total(indic), by(year)
	egen denom = total(not), by(year)
	gen stat = (numer/denom)*100

	preserve

	duplicates drop year stat, force
	sort year

	twoway line stat year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
	title(Percentage of Repeat Sales by Property) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_cre_repeat_prop.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_repeat_prop.dta, replace

	restore


	****Here we collect the percentage change in prices for repeat sales by property
	duplicates drop g_pkid Status_dt, force


	****Create the price change variable
	gen price_change = .
	bysort g_pkid (Status_dt): replace price_change = (Price[_n]/Price[_n-1]) - 1


	****Find the price changes for the repeat sales, their annual means, and graph both statistics
	gen pc_repeat = (price_change*indic)*100
	egen pcr_mean = mean(pc_repeat), by(year)


	****Scatter
	format Status_dt %tdmonCCYY
	twoway scatter pc_repeat Status_dt, mcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Date) ytitle(Percent) ///
	title(Change in Price for Repeat Sales) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_repeat_scatter.png, replace

	saveold ~/analysts/tanya/`poi'_repeat_scatter.dta, replace
	
	format Status_dt %td

	preserve
	duplicates drop pcr_mean year, force
	sort year


	****Time Series of annual means
	twoway line pcr_mean year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
	title(Average Change in Price for Repeat Sales) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_cre_repeat_ch_pr_mean.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_repeat_ch_pr_mean.dta, replace

	restore


	****Now we look at the breakdown of sales by property type
	levelsof MainType, local(typez)
	local tic = 1


	****Create percentage variables for each property type
	foreach x in `typez' {
		gen type_`tic' = 0
		replace type_`tic' = 1 if MainType == "`x'"
		egen ty_to_`tic' = total(type_`tic'), by(year)
		gen type_stat_`tic' = (ty_to_`tic'/kp)*100
		local tic = `tic' + 1
		}


	****Labels for the graphs
	local dtype_stat_1 "Apartment"
	local dtype_stat_2 "Dev Site"
	local dtype_stat_3 "Hotel"
	local dtype_stat_4 "Industrial"
	local dtype_stat_5 "Office"
	local dtype_stat_6 "Retail"


	****One graph for each type
	forvalues x = 1/6 {
		preserve
		duplicates drop type_stat_`x' year, force
		sort year
		twoway line type_stat_`x' year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
		title("Percent of Sales by Property Type: `dtype_stat_`x''") subtitle("`d`poi''")
		graph export ~/analysts/tanya/`poi'_type_stat_`x'.png, replace
		saveold  ~/analysts/tanya/`poi'_type_stat_`x'.dta, replace
		restore
		}
		

	****Now we find out how many of the new sales are potentially troubled and distress assets
	merge m:1 PropertyKey_id using ~/analysts/tanya/ptdassets1013.dta

	gen ptd = 0
	replace ptd = 1 if _m == 3
	drop _m

	levelsof g_pkid if ptd == 1, local(ptd_pkid)
	drop indic
	gen indic = .

	foreach b in `ptd_pkid' {	
		levelsof Status_dt if g_pkid == `b', local(p_times)
		foreach y in `p_times' {
			local plus = 0
			sum g_pkid if (Troubled_dt >= `y' & g_pkid == `b')
			local plus = `plus' + `r(N)'
			replace indic = 1 if (g_pkid == `b' & Status_dt == `y' & `plus' != 0)
			}
		}

	replace indic = 0 if indic >= .


	****Now we graph potentially troubled and distressed assets as a percentage of new sales
	drop not numer denom stat
	gen not = 1
	egen numer = total(indic), by(year)
	egen denom = total(not), by(year)
	gen stat = (numer/denom)*100

	preserve

	duplicates drop year stat, force
	sort year

	twoway line stat year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
	title(PTD Assets as a Percentage of Sales) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_cre_ptd_prop.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_ptd_prop.dta, replace

	restore


	****Now we graph dollar amounts of potentially troubled and distressed assets as a percentage of new sales
	drop not numer denom stat
	gen not = 1

	egen numer = total(indic*Price), by(year)
	egen denom = total(not*Price), by(year)
	gen stat = (numer/denom)*100

	preserve

	duplicates drop year stat, force
	sort year

	twoway line stat year, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
	title(PTD Assets as a Percentage of Sales: Dollar Amounts) subtitle("`d`poi''")
	graph export ~/analysts/tanya/`poi'_cre_ptd_prop_dollar.png, replace
	
	saveold ~/analysts/tanya/`poi'_cre_ptd_prop_dollar.dta, replace

	restore

	}
	

****Here we graph the unsold REO PTD assets as a percentage of PTD assets
use ~/analysts/tanya/FR1013.dta, clear

merge m:1 PropertyKey_id using ~/analysts/tanya/ptdassets1013.dta

gen ptd = 0
replace ptd = 1 if _m == 3
drop _m

gen kp = 0
replace kp = 1 if ((ptd == 1 & Troubled_dt > Status_dt) | (ptd == 0))
keep if kp == 1

gen indic = 1
gen year_ptd = year(Troubled_dt)

gen reo = 0
replace reo = 1 if DistressedStatus_tx == "LenderREO"

egen numer = total(reo), by(year_ptd)
egen denom = total(indic), by(year_ptd)

gen stat = (numer/denom)*100

preserve

duplicates drop year_ptd stat, force
sort year_ptd

twoway line stat year_ptd, lwidth(thick) lpattern(dash_dot) lcolor(blue) graphregion(color(white)) ylabel(, nogrid) xtitle(Year) ytitle(Percent) ///
title(REO PTD Assets as a Percentage of PTD Assets) subtitle(Unsold Properties)
graph export ~/analysts/tanya/cre_ptd_reo.png, replace

saveold ~/analysts/tanya/cre_ptd_reo.dta, replace

restore




