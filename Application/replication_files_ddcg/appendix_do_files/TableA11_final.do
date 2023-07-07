/**********************/
/*Clean and set memory*/
/**********************/
clear all
set matsize 5000
set more off

/******************/
/*Install programs*/
/******************/

/*********************/
/*Sets base directory*/
/*********************/
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */

/***************************************/
/*Analysis of propensity to democratize*/
/***************************************/
use ${project}/DDCGdata_final.dta, clear
replace y=y/100
gen dy=d.y

gen tdemoc=.
replace tdemoc=1 if dem==1&l.dem==0
replace tdemoc=0 if dem==0&l.dem==0

/* Estimates */
probit tdemoc l.y i.year, r cluster(wbcode2)
margins,  dydx(L1.y) post
nlcom (inceffect: _b[L1.y]), post
estimates store e1

probit tdemoc l.dy l2.y yy*, r
margins,  dydx(L1.dy L2.y) post
nlcom (dif1: _b[L1.dy]) (inceffect: _b[L2.y]), post
estimates store e2

probit tdemoc l.dy l2.dy l3.y yy*, r
margins,  dydx(L1.dy L2.dy L3.y) post
nlcom (dif1: _b[L1.dy]) (dif2: _b[L2.dy]) (inceffect: _b[L3.y]), post
estimates store e3

probit tdemoc l.dy l2.dy l3.dy l4.y yy*, r
predict _pscore1 if e(sample), pr 
margins,  dydx(L1.dy L2.dy L3.dy L4.y) post
nlcom (dif1: _b[L1.dy]) (dif2: _b[L2.dy]) (dif3: _b[L3.dy]) (inceffect: _b[L4.y]), post
estimates store e4

probit tdemoc l.dy l2.dy l3.dy  yy*, r
predict _pscore2 if e(sample), pr 
margins,  dydx(L1.dy L2.dy L3.dy) post
nlcom (dif1: _b[L1.dy]) (dif2: _b[L2.dy]) (dif3: _b[L3.dy]) , post
estimates store e5

/*Analysis of pscore*/
sum _pscore1 _pscore2
cor _pscore1 _pscore2

/*Overlap figures*/
kdensity _pscore1 if tdemoc==1, generate(x_treat1 prob_treat1) nograph
kdensity _pscore1 if tdemoc==0, generate(x_control1 prob_control1) nograph

twoway (line prob_treat1 x_treat1, lcolor(black) lpattern(solid)) ///
       (line prob_control1 x_control1, lcolor(gray) lpattern(dash)) ///
       , ytitle(Estimated density) xtitle(Propensity to democratize) legend(on order(1 "Propensity to democratize among countries that democratized" 2 "Propensity to democratize among countries that remained nondemocratic" ) rows(2) region(fcolor(white) margin(zero) lcolor(white)) bmargin(zero)) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/OverlapDemoc1.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/

/*Exports tables*/
estout e1 e2 e3 e4 e5  using ${project}/appendix/TablePscoreTop.tex, ///
style(tex) varlabels(dif1 " Change in GDP at $ t-1$ " dif2 " Change in GDP at $ t-2$ " dif3 " Change in GDP at $ t-3$ " inceffect " GDP level effect ") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N N_g, fmt(%7.0f %7.0f) labels( " Observations" "Countries in sample"))   collabels(none) ///
keep(dif1 dif2 dif3 inceffect) order(dif1 dif2 dif3 inceffect) stardrop(dif1 dif2 dif3 inceffect) nolabel replace mlabels(none)


/*************************************************/
/*Analysis of propensity to revert from democracy*/
/*************************************************/
use ${project}/DDCGdata_final.dta, clear
replace y=y/100
gen dy=d.y

gen trever=.
replace trever=1 if dem==0&l.dem==1
replace trever=0 if dem==1&l.dem==1

/* Estimates */
probit trever l.y yy*, r
margins,  dydx(L1.y) post
nlcom (inceffect: _b[L1.y]), post
estimates store e1b

probit trever l.dy l2.y yy*, r
margins,  dydx(L1.dy L2.y) post
nlcom (dif1: _b[L1.dy]) (inceffect: _b[L2.y]), post
estimates store e2b

probit trever l.dy l2.dy l3.y yy*, r
margins,  dydx(L1.dy L2.dy L3.y) post
nlcom (dif1: _b[L1.dy]) (dif2: _b[L2.dy]) (inceffect: _b[L3.y]), post
estimates store e3b

probit trever l.dy l2.dy l3.dy l4.y yy*, r
predict _pscore1 if e(sample), pr 
margins,  dydx(L1.dy L2.dy L3.dy L4.y) post
nlcom (dif1: _b[L1.dy]) (dif2: _b[L2.dy]) (dif3: _b[L3.dy]) (inceffect: _b[L4.y]), post
estimates store e4b

probit trever l.dy l2.dy l3.dy  yy*, r
predict _pscore2 if e(sample), pr 
margins,  dydx(L1.dy L2.dy L3.dy) post
nlcom (dif1: _b[L1.dy]) (dif2: _b[L2.dy]) (dif3: _b[L3.dy]), post
estimates store e5b


/*Analysis of pscore*/
sum _pscore1 _pscore2
cor _pscore1 _pscore2

/*Overlap figures*/
kdensity _pscore1 if trever==1, generate(x_treat1 prob_treat1) nograph
kdensity _pscore1 if trever==0, generate(x_control1 prob_control1) nograph

twoway (line prob_treat1 x_treat1, lcolor(black) lpattern(solid)) ///
       (line prob_control1 x_control1, lcolor(gray) lpattern(dash)) ///
       , ytitle(Estimated density) xtitle(Propensity to revert) legend(on order(1 "Propensity to revert among countries that reverted" 2 "Propensity to revert among countries that remained in democracy" ) rows(2) region(fcolor(white) margin(zero) lcolor(white)) bmargin(zero)) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/OverlapRever1.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/

/*Exports tables*/
estout e1b e2b e3b e4b e5b  using ${project}/appendix/TablePscoreBot.tex, ///
style(tex) varlabels(dif1 " Change in GDP at $ t-1$ " dif2 " Change in GDP at $ t-2$ " dif3 " Change in GDP at $ t-3$ " inceffect " GDP level effect ") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N N_g, fmt(%7.0f %7.0f) labels( " Observations" "Countries in sample"))   collabels(none) ///
keep(dif1 dif2 dif3 inceffect) order(dif1 dif2 dif3 inceffect) stardrop(dif1 dif2 dif3 inceffect) nolabel replace mlabels(none)
