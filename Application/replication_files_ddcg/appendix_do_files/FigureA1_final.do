/**********************/
/*Clean and set memory*/
/**********************/
clear all
set matsize 5000
set more off

/******************/
/*Install programs*/
/******************/
cap ssc install xtabond2 
cap ssc install xtivreg2 
cap ssc install spmat 
cap ssc install spmack

/*********************/
/*Sets base directory*/
/*********************/
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */

/****************************/
/***Democracy indices****/
/****************************/

use ${project}/DDCGdata_final.dta, clear

gen pol4=polity if polity>=-10
replace pol4=(pol4+10)/20
gen fh=(14-cl-pr)/12

/*Figures by regions*/

preserve 
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in the World) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemWorld.eps, as(eps) preview(off) replace
restore



preserve 
keep if region=="AFR"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in Africa) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemAFR.eps, as(eps) preview(off) replace
restore

preserve 
keep if region=="EAP"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
 ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in East Asia and the Pacific) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemEAP.eps, as(eps) preview(off) replace
restore

preserve 
keep if region=="ECA"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
 ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in Eastern Europe and Central Asia) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemECA.eps, as(eps) preview(off) replace
restore

preserve 
keep if region=="INL"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
 ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in Western Europe and Offshoots) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemINL.eps, as(eps) preview(off) replace
restore

preserve 
keep if region=="LAC"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
 ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in Latin American and the Caribbean) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemLAC.eps, as(eps) preview(off) replace
restore

preserve 
keep if region=="MNA"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"
label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
 ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in Middle East and North Africa) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemMNA.eps, as(eps) preview(off) replace
restore

preserve 
keep if region=="SAS"
collapse dem demPS demBMR demCGV pol4 fh, by(year)
tsset year
label variable dem "Our democracy measure"
label variable demPS "PS democracy measure"
label variable demCGV "CGV democracy measure"
label variable demBMR "BMR democracy measure"

label variable pol4 "Polity IV index"
label variable fh "Freedom House index"

twoway (tsline dem, lcolor(red)) ///
(tsline demPS, lcolor(black) lpattern(dash)) ///
(tsline demCGV, lcolor(black) lpattern(dot)) ///
(tsline demBMR, lcolor(black) lpattern(longdash_dot)) ///
(tsline pol4, lcolor(blue) lpattern(solid)) ///
(tsline fh, lcolor(blue) lpattern(dash)), ///
 ytitle(Democracy index average) xtitle(Year) title(Democracy 1960-2010 in South Asia) xlabel(1960(10)2010) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export ${project}/appendix/FigureDemSAS.eps, as(eps) preview(off) replace
restore
