/**********************/
/*Clean and set memory*/
/**********************/
clear all
set matsize 5000
set more off

/******************/
/*Install programs*/
/******************/
cap ssc install parmest

/*********************/
/*Sets base directory*/
/*********************/
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/*********************DEFINE REQUIRED PROGRAMS THAT WILL BE USED DURING THE EXECUTION OF THIS DO FILE **************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
program tefpas_ra, eclass

quietly: tab year if free_elections_Lag15!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (free_elections_Lag15 lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=_b[ATET: r1vs0.tdemoc]
drop  dumyears*

quietly: tab year if constrained_Lag15!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (constrained_Lag15 lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if inclusive_Lag15!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (inclusive_Lag15 lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if civil_rights_Lag15!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (civil_rights_Lag15 lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if dem_Lag15!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (dem_Lag15 lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

local column_names  free_elections_time0 constrained_time0 inclusive_time0 civil_rights_time0 dem_time0

forvalues s=14(-1)2{
display `s'

quietly: tab year if free_elections_Lag`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (free_elections_Lag`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if constrained_Lag`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (constrained_Lag`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if inclusive_Lag`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (inclusive_Lag`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if civil_rights_Lag`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (civil_rights_Lag`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if dem_Lag`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (dem_Lag`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

local num=15-`s'
local column_names `column_names' free_elections_time`num' constrained_time`num' inclusive_time`num' civil_rights_time`num' dem_time`num'
}

matrix b=(b, 0, 0, 0, 0, 0)
local column_names `column_names' free_elections_time15 constrained_time15 inclusive_time15 civil_rights_time15 dem_time15

forvalues s=0(1)20{
display `s'

quietly: tab year if free_elections_Lead`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (free_elections_Lead`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if constrained_Lead`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (constrained_Lead`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if inclusive_Lead`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (inclusive_Lead`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if civil_rights_Lead`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (civil_rights_Lead`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

quietly: tab year if dem_Lead`s'!=. & tdemoc!=., gen(dumyears)
quietly: teffects ra (dem_Lead`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (tdemoc), atet iterate(5)
matrix b=(b, _b[ATET: r1vs0.tdemoc])
drop  dumyears*

local num=16+`s'
local column_names `column_names' free_elections_time`num' constrained_time`num' inclusive_time`num' civil_rights_time`num' dem_time`num' 
}

mat colnames b = `column_names'
ereturn post b
end


/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/************************************ ESTIMATION PROCEDURES AND STORING RESULTS ************************************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/

use ${project}/DDCGdata_final.dta

sort wbcode2 year

gen tdemoc=.
replace tdemoc=1 if dem==1&l.dem==0
replace tdemoc=0 if dem==0&l.dem==0

gen lag1y=l1.y
gen lag2y=l2.y
gen lag3y=l3.y
gen lag4y=l4.y

keep if lag1y!=.
keep if lag2y!=.
keep if lag3y!=.
keep if lag4y!=.

tab year, gen(dyear)


gen civil_rights=(7-cl)/6

gen free_elections=0 if xrcomp!=.
replace free_elections=1 if (xrcomp==2|xrcomp==3)&(xropen==3|xropen==4)

gen inclusive=0 if parcom!=.
replace inclusive=1 if parcomp==4|parcomp==5

gen constrained=0 if xconst!=.
replace constrained=1 if xconst==5|xconst==6|xconst==7


forvalues j=15(-1)1{
gen  free_elections_Lag`j'=L`j'.free_elections-L.free_elections
gen  constrained_Lag`j'=L`j'.constrained-L.constrained
gen  inclusive_Lag`j'=L`j'.inclusive-L.inclusive
gen  civil_rights_Lag`j'=L`j'.civil_rights-L.civil_rights
gen  dem_Lag`j'=L`j'.dem-L.dem
}

forvalues j=0(1)20{
gen  free_elections_Lead`j'=F`j'.free_elections-L.free_elections
gen  constrained_Lead`j'=F`j'.constrained-L.constrained
gen  inclusive_Lead`j'=F`j'.inclusive-L.inclusive
gen  civil_rights_Lead`j'=F`j'.civil_rights-L.civil_rights
gen  dem_Lead`j'=F`j'.dem-L.dem
}

keep if tdemoc!=.
set seed 12345

xtset, clear
bootstrap _b, reps(50) cluster(wbcode2): tefpas_ra 
parmest, format(estimate min95 max95) saving(${project}/appendix/impulse_components, replace)

/*Plots: elections*/
use ${project}/appendix/impulse_components, clear
split parm, p("free_elections_time")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Free elections dummy") xtitle(Years around democratization) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)20) 
graph export ${project}/appendix/components_FreeElections.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/

/*Plots: constraints*/
use ${project}/appendix/impulse_components, clear
split parm, p("constrained_time")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Constraints on the executive dummy") xtitle(Years around democratization) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)20) 
graph export ${project}/appendix/components_Constraints.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/

/*Plots: inclusive*/
use ${project}/appendix/impulse_components, clear
split parm, p("inclusive_time")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Inclussive politics dummy") xtitle(Years around democratization) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)20) 
graph export ${project}/appendix/components_Inclussive.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/


/*Plots: civil liberties*/
use ${project}/appendix/impulse_components, clear
split parm, p("civil_rights_time")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Civil liberties") xtitle(Years around democratization) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)20) 
graph export ${project}/appendix/components_Civil.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/

/*Plots: democracy*/
use ${project}/appendix/impulse_components, clear
split parm, p("dem_time")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Democracy persistence") xtitle(Years around democratization) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)20) 
graph export ${project}/appendix/components_Democracy.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/






