/* Project:	CVA
   Input Path:	XXXX/
   Output Path:	/data8/home/l1mdr03/CVA/	
   Date:	19 September 2012
*/

/*The data used in the following program was prepared by Hamed Faquiryan and is a merge of the Y9-C, Compustat, and
Markit CDS datasets using the CRSP/Y9-C map from the NY Fed, the CRSP/Compustat link-table from the BOG,
and hand-matching CDS tickers. The program seeks to estimate the effects of CVA adjustments on trading behavior.
The analysis conducted here is for Lopez & Rodriguez (2013).*/

clear all
set more off
cd XXXX/CVA
log using cva_log, text replace

use XXXX/CVA/cva_10_5.dta


/*To begin, we first generate log-growth in trading revenue and drop any missing observations,
which is a large subset of the data. Then we keep only the firms with more than 10 Billion assets*/

sort gvkey date
by gvkey: gen grev=log(trade_rev[_n]/trade_rev[_n-1])
drop if grev==.
by gvkey: egen massets=mean(tot_assets)
drop if massets<10000000
drop if gvkey==3555 & scfq==. /*Countrywide has a repeat observation for 2002q1*/
xtset gvkey date, q
save bhc_derv_vars, replace


/*Now we'll create a rolling volatility measure for grev with a 1.5 year window and merge with the above.*/
 
quietly rolling sdrev=r(sd), window(6) keep(gvkey) clear: summarize grev
drop date start
rename end date
sort gvkey date
merge 1:1 gvkey date using bhc_derv_vars
keep if _merge==3
drop if sdrev==.
sort gvkey date
drop _merge
erase bhc_derv_vars.dta


//We now run a fixed effects on growth of revenues to create a volatility estimate
quietly xi: reg grev i.gvkey i.date, robust
predict grevh
gen res=grev-grevh
replace res=abs(res)
label var res "RESVOL"
label var sdrev "SDREV"
drop _Igvkey* _Idate*
**keep if bhckk090!=.	/*Uncomment this piece if you want only firms with adjustment data*/


/*Now we generate a series of control variables*/

gen totpfv=fv_pirderivs+fv_pcred+fv_pfxderivs+fv_peqderivs+fv_pcomderivs
label var totpfv "Total Positive Fair Value"

gen totnfv=fv_nirderivs+fv_ncred+fv_nfxderivs+fv_neqderivs+fv_ncomderivs
label var totnfv "Total Negative Fair Value"

gen dacon=(fv_pirderivs^2+fv_pcred^2+fv_pfxderivs^2+fv_peqderivs^2+fv_pcomderivs^2)/totpfv^2
label var dacon "Derivative Asset Concentration"

gen dlcon=(fv_nirderivs^2+fv_ncred^2+fv_nfxderivs^2+fv_neqderivs^2+fv_ncomderivs^2)/totnfv^2
label var dlcon "Derivative Liability Concentration"

gen fv_irp=fv_pirderivs/totpfv
label var fv_irp "Share of IR Derivative Asset"

gen fv_credp=fv_pcred/totpfv
label var fv_credp "Share of Credit Derivative Asset"

gen fv_fxp=fv_pfxderivs/totpfv
label var fv_fxp "Share of FX Derivative Asset"

gen fv_eqp=fv_peqderivs/totpfv
label var fv_eqp "Share of Equities Derivative Asset"

gen fv_comp=fv_pcomderivs/totpfv
label var fv_comp "Share of Commodities Derivative Asset"

gen fv_irn=fv_nirderivs/totnfv
label var fv_irn "Share of IR Derivative Liability"

gen fv_credn=fv_ncred/totnfv
label var fv_credn "Share of Credit Derivative Liability"

gen fv_fxn=fv_nfxderivs/totnfv
label var fv_fxn "Share of FX Derivative Liability"

gen fv_eqn=fv_neqderivs/totnfv
label var fv_eqn "Share of Equities Derivative Liability"

gen fv_comn=fv_ncomderivs/totnfv
label var fv_comn "Share of Commodities Derivative Liability"

gen lirp=log(1+fv_pirderivs)
gen lcredp=log(1+fv_pcred)
gen lfxp=log(1+fv_pfxderivs)
gen leqp=log(1+fv_peqderivs)
gen lcomp=log(1+fv_pcomderivs)
gen lirn=log(1+fv_nirderivs)
gen lcredn=log(1+fv_ncred)
gen lfxn=log(1+fv_nfxderivs)
gen leqn=log(1+fv_neqderivs)
gen lcomn=log(1+fv_ncomderivs)
gen lbanks=log(1+ncrd_expsr_banks)
gen lfinguars=log(1+ncrd_expsr_finguars)
gen lhdgfnds=log(1+ncrd_expsr_hdgfnds)
gen lsovgovs=log(1+ncrd_expsr_sovgovs)
gen lcorps=log(1+ncrd_expsr_corps)

gen VaR=log(1+VaRf/12.5)
label var VaR "VaR"

gen adjsr=(cva+dva)/trade_rev
label var adjsr "Total Adjustment to Trading Revenue"

gen nadjsrevr=1-adjsr
gen firevr=ir_expsr/trade_rev
gen fxrevr=fx_expsr/trade_rev
gen eqrevr=eq_expsr/trade_rev
gen comrevr=com_expsr/trade_rev
gen crrevr=crd_expsr_rev/trade_rev

gen size=log(1+trade_assets)
label var size "Trading Assets to Total Assets"

gen cvar=cva/trade_rev
label var cvar "Counterparty Adjustment to Trading Revenue"

gen dvar=dva/trade_rev
label var dvar "Own Adjustment to Trading Revenue"

gen tot_ncrd=ncrd_expsr_banks+ncrd_expsr_finguars+ncrd_expsr_hdgfnds+ncrd_expsr_sovgovs+ncrd_expsr_corps
gen tot_col=totcol_banks+totcol_finguars+totcol_hdgfnds+totcol_sovgovs+totcol_corps

gen banksr=totcol_banks/ncrd_expsr_banks
label var banksr "Total Collateral to Net Exposure - Banks"

gen finguarsr=totcol_finguars/ncrd_expsr_finguars
label var finguarsr "Total Collateral to Net Exposure - Financial Guarantors"

gen hdgfndsr=totcol_hdgfnds/ncrd_expsr_hdgfnds
label var hdgfndsr "Total Collateral to Net Exposure - Hedge Funds"

gen sovgovsr=totcol_sovgovs/ncrd_expsr_sovgovs
label var sovgovsr "Total Collateral to Net Exposure - Sovereigns"

gen corpsr=totcol_corps/ncrd_expsr_corps
label var corpsr "Total Collateral to Net Exposure - Corporates"

gen colr=tot_col/tot_ncrd
label var colr "Ratio of Total Collateral to Total Net Current Exposure"

label var spread5y30mn "Mean of Log 5Y CDS Spread"

save XXXX/CVA/y9c_adjustments, replace	

 
**Now we'll create the Slope and default variables that will be used in the regressions
drop _all
insheet using  XXXX/Swap_Daily.txt, t
gen da=date(date,"YMD")
gen month=month(da)
gen year=year(da)
gen day=day(da)
gen quarter=0
replace quarter=1 if month==1|month==2|month==3
replace quarter=2 if month==4|month==5|month==6
replace quarter=3 if month==7|month==8|month==9
replace quarter=4 if month==10|month==11|month==12
gen yq=yq(year,quarter)
gen slp=log(dswp10/dswp1)
drop if dswp1==.
sort yq month day
by yq: keep if _n==_N
gen cslp=slp[_n]-slp[_n-1]
keep yq cslp
drop if cslp==.
rename yq date
sort date
label var cslp "Change in Swap Curve Slope"
save XXXX/CVA/Slope, replace


**Now the VIX data
drop _all
insheet using XXXX/Rates_Daily_Close.txt, t
rename vixcls vix
gen da=date(date,"YMD")
gen month=month(da)
gen year=year(da)
gen day=day(da)
gen quarter=0
replace quarter=1 if month==1|month==2|month==3
replace quarter=2 if month==4|month==5|month==6
replace quarter=3 if month==7|month==8|month==9
replace quarter=4 if month==10|month==11|month==12
gen yq=yq(year,quarter)
sort yq day
by yq: keep if _n==_N
gen clvix=log(vix[_n]/vix[_n-1])
gen cvix=(vix[_n]-vix[_n-1])/vix[_n-1]
gen lvix=log(vix)
keep yq lvix clvix cvix vix
rename yq date
merge 1:1 date using Slope
keep if _merge==3
drop _merge
sort date
label var lvix "Log of VIX"
label var clvix "Log Change of VIX"
save XXXX/CVA/Slope, replace


**Now we'll merge the financial variables datasets with the bank dataset
drop _all
use XXXX/CVA/y9c_adjustments
sort date
merge m:1 date using Slope
keep if _merge==3|_merge==1
drop _merge
sort entity date
save XXXX/CVA/y9c_cva, replace


**Finally, creating the CDS indices data
drop _all
use XXXX/cds_indices_7_20.dta
replace invrgd_cdx=invrgd_cdx/10000
replace hyld_cdx=hyld_cdx/10000
replace sov_cdx=sov_cdx/10000
gen year=year(date)
gen month=month(date)
gen q=.
replace q=1 if month==1|month==2|month==3
replace q=2 if month==4|month==5|month==6
replace q=3 if month==7|month==8|month==9
replace q=4 if month==10|month==11|month==12
gen ln_invgr=ln(invrgd_cdx)
gen ln_hiyld=ln(hyld_cdx)
gen ln_sov=ln(sov_cdx)
gen yq=yq(year,q)
sort yq date
by yq: keep if _n==_N
sort yq
drop date
rename yq date
save cdsind, replace
drop _all
use XXXX/CVA/y9c_cva
sort date
merge m:1 date using cdsind
keep if _merge==3|_merge==1
drop _merge
sort gvkey date
save XXXX/CVA/y9c_cva, replace


**Last-minute clean-up

gen anoqr=anoq/atq
gen lnoqr=lnoq/ltq


**Summary Statitics for the variables included in the regressions
estpost summarize sdrev res adjsr cvar dvar size dacon dlcon VaR ///
fv_irn fv_credn fv_fxn fv_eqn fv_comn fv_irp fv_credp fv_fxp fv_eqp fv_comp, detail
esttab . using sum_y9c.tex, cells("mean(fmt(4)) p50(fmt(4)) sd(fmt(4))") noobs nonumber booktabs ///
wide title(Summary Statistics\label{sumy9c}) replace


**Correlation matrix
estpost correlate fv_irn fv_credn fv_fxn fv_eqn fv_comn fv_irp fv_credp fv_fxp fv_eqp fv_comp ///
spread5y30mn ln_invgr ln_hiyld, matrix
esttab . using correlation2.tex, unstack b(4) nonotes noobs nostar booktabs compress ///
title(Correlation Matrix) nonumber replace 

**Second correlation
estpost correlate sdrev res adjsr cvar dvar size dacon dlcon VaR, matrix
esttab . using correlation1.tex, unstack b(4) nonotes noobs nostar booktabs compress ///
title(Correlation Matrix) nonumber replace 


/*Now the regression results*/

**The first set of regressions look at the relationship between Trading Revenue Volatility and adjustments
eststo clear
eststo: quietly xtreg sdrev adjsr, fe robust
eststo: quietly xtreg sdrev adjsr size, fe robust
eststo: quietly xtreg sdrev adjsr size dacon dlcon, fe robust
eststo: quietly xtreg sdrev adjsr size dacon dlcon VaR, fe robust
eststo: quietly xtreg sdrev adjsr size dacon dlcon VaR cslp clvix, fe robust
eststo: quietly xtreg sdrev size dacon dlcon VaR cslp clvix, fe robust
esttab using sdrev_reg.tex, se ar2 scalars(N_g F) booktabs title(SDREV Regression\label{sdrevreg}) replace


/*The second set of regressions look at the relationship between Trading Revenue Residual Volatility and adjustments*/

eststo clear
eststo: quietly xtreg res adjsr, fe robust
eststo: quietly xtreg res adjsr size, fe robust
eststo: quietly xtreg res adjsr size dacon dlcon, fe robust
eststo: quietly xtreg res adjsr size dacon dlcon VaR, fe robust
eststo: quietly xtreg res adjsr size dacon dlcon VaR cslp clvix, fe robust
eststo: quietly xtreg res size dacon dlcon VaR cslp clvix, fe robust
esttab using res_reg.tex, se ar2 scalars(N_g F) booktabs title(RES Regression\label{resreg}) replace


/*The third Set of Regressions look at the relationship between Trading Revenue Volatility and detailed adjustments*/
eststo clear
eststo: quietly xtreg sdrev cvar dvar, fe robust
eststo: quietly xtreg sdrev cvar dvar size, fe robust
eststo: quietly xtreg sdrev cvar dvar size dacon dlcon, fe robust
eststo: quietly xtreg sdrev cvar dvar size dacon dlcon VaR, fe robust
eststo: quietly xtreg sdrev cvar dvar size dacon dlcon VaR cslp clvix, fe robust
eststo: quietly xtreg sdrev size dacon dlcon VaR cslp clvix, fe robust
esttab using sdrev_regd.tex, se ar2 scalars(N_g F) booktabs title(SDREV Detailed Regression\label{sdrevregd}) ///
replace


/*The fourth set of Regressions look at the relationship between Trading Revenue Residual Volatility and detailed adjustments*/

eststo clear
eststo: quietly xtreg res cvar dvar, fe robust
eststo: quietly xtreg res cvar dvar size, fe robust
eststo: quietly xtreg res cvar dvar size dacon dlcon, fe robust
eststo: quietly xtreg res cvar dvar size dacon dlcon VaR, fe robust
eststo: quietly xtreg res cvar dvar size dacon dlcon VaR cslp clvix, fe robust
eststo: quietly xtreg res size dacon dlcon VaR cslp clvix, fe robust
esttab using res_regd.tex, se ar2 scalars(N_g F) booktabs title(RES Detailed Regression\label{resregd}) replace


/*The fifth regression accounts for the DVA adjustments - Y9-C. Also, we need some extra controls.*/

eststo clear
eststo: quietly xtreg dvar spread5y30mn, fe robust
eststo: quietly xtreg dvar spread5y30mn fv_irn fv_credn fv_fxn fv_eqn , fe robust
eststo: quietly xtreg dvar spread5y30mn fv_irn fv_credn fv_fxn fv_eqn  dlcon, fe robust
esttab using dvar_reg.tex, se ar2 scalars(N_g F) booktabs title(DVAR Regression\label{dvarreg}) replace


/*The sixth regression accounts for the CVA adjustments - Y9-C*/

eststo clear
eststo: quietly xtreg cvar ln_invgr, fe robust
eststo: quietly xtreg cvar ln_invgr fv_irp fv_credp fv_fxp fv_eqp  dacon, fe robust
eststo: quietly xtreg cvar ln_hiyld, fe robust
eststo: quietly xtreg cvar ln_hiyld fv_irp fv_credp fv_fxp fv_eqp  dacon, fe robust
esttab using cvar_reg.tex, se ar2 scalars(N_g F) booktabs title(CVAR Regression\label{cvarreg}) replace


**Now for some beautiful graphs
xtset gvkey date, q
drop year
gen year = year(dofq(date))

**XXXX Volatilities
twoway (line sdrev date if gvkey==2968&year>=2007 , yaxis(1)) ///
||(line vix date if gvkey==2968&year>=2007, yaxis(2)),title("XXXX") ///
ytitle("VIX",axis(2)) ytitle("Vol. Measure - SDREV") scheme(s2mono) ///
legend(label(1 "Std. Dev. Trad. Rev.") label(2 "VIX")) graphregion(color(white))
graph export jpmc_sdrev_plot.png, replace
twoway (line res date if gvkey==2968&year>=2007 , yaxis(1)) ///
||(line vix date if gvkey==2968&year>=2007, yaxis(2)),title("XXXX") ///
ytitle("VIX",axis(2)) ytitle("Vol. Measure - RES") scheme(s2mono) ///
legend(label(1 "Vol. Trad. Rev.") label(2 "VIX")) graphregion(color(white))
graph export jpmc_res_plot.png, replace

**XXXX Volatilities
twoway (line sdrev date if gvkey==7647&year>=2007 , yaxis(1)) ///
||(line vix date if gvkey==7647&year>=2007, yaxis(2)),title("XXXX") ///
ytitle("VIX",axis(2)) ytitle("Vol. Measure - SDREV") scheme(s2mono) ///
legend(label(1 "Std. Dev. Trad. Rev.") label(2 "VIX")) graphregion(color(white))
graph export bac_sdrev_plot.png, replace
twoway (line res date if gvkey==7647&year>=2007 , yaxis(1)) ///
||(line vix date if gvkey==7647&year>=2007, yaxis(2)),title("XXXX") ///
ytitle("VIX",axis(2)) ytitle("Vol. Measure - RES") scheme(s2mono) ///
legend(label(1 "Vol. Trad. Rev.") label(2 "VIX")) graphregion(color(white))
graph export bac_res_plot.png, replace


**Now to adjustment graphs
gen cvag=cva/1000000
gen dvag=dva/1000000

**XXXX Adjustments
twoway (line cvag dvag date if gvkey==2968&year>=2011, yaxis(1)) ///
||(line spread5y30mn date if gvkey==2968&year>=2011, yaxis(2)),title("XXXX Impact on Trading Revenue ") ///
ytitle("CDS 5Y Spread",axis(2)) ytitle("Billions of US$") scheme(s2mono) ///
legend(label(1 "CVA") label(2 "DVA") label(3 "CDS")) graphregion(color(white))
graph export jpmc_cva_plot.png, replace

**XXXX Adjustmets
twoway (line cvag dvag date if gvkey==7647&year>=2011, yaxis(1)) ///
||(line spread5y30mn date if gvkey==7647&year>=2011, yaxis(2)),title("XXXX Impact on Trading Revenue ") ///
ytitle("CDS 5Y Spread",axis(2)) ytitle("Billions of US$") scheme(s2mono) ///
legend(label(1 "CVA") label(2 "DVA") label(3 "CDS")) graphregion(color(white))
graph export bac_cva_plot.png, replace

**Scatter Plot of DVA and CVA
twoway (scatter cvag dvag)||(lfitci cvag dvag), ytitle("CVA Billions of US$") ///
xtitle("DVA Billions of US$") title("CVA vs DVA") ///
scheme(s2mono) graphregion(color(white))
graph export cva_dva_plot.png, replace


log close


/*Aaaaaaaaannnnddddd we're done!*/

