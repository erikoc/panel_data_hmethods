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
global limit=25                                                                             /* Evaluate effects 25 years after transition */
mata: mata set matafavor speed

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/*********************DEFINE REQUIRED PROGRAMS THAT WILL BE USED DURING THE EXECUTION OF THIS DO FILE **************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/

capture program drop vareffects
program define vareffects, eclass

quietly: nlcom (effect1: _b[shortrun]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post

quietly: nlcom (effect2: _b[effect1]*_b[lag1]+_b[shortrun]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post

quietly: nlcom (effect3: _b[effect2]*_b[lag1]+_b[effect1]*_b[lag2]+_b[shortrun]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
	  
quietly: nlcom (effect4: _b[effect3]*_b[lag1]+_b[effect2]*_b[lag2]+_b[effect1]*_b[lag3]+_b[shortrun]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post	  

forvalues j=5(1)$limit{	  
local j1=`j'-1
local j2=`j'-2
local j3=`j'-3
local j4=`j'-4

quietly: nlcom (effect`j': _b[effect`j1']*_b[lag1]+_b[effect`j2']*_b[lag2]+_b[effect`j3']*_b[lag3]+_b[effect`j4']*_b[lag4]+_b[shortrun]) ///
	  (effect`j1': _b[effect`j1']) ///
	  (effect`j2': _b[effect`j2']) ///
	  (effect`j3': _b[effect`j3']) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post	  	  

}

quietly: nlcom (effect$limit: _b[effect$limit]) ///
	  (longrun: _b[shortrun]/(1-_b[lag1]-_b[lag2]-_b[lag3]-_b[lag4])) ///
      (shortrun: _b[shortrun]) ///
	  (persistence: _b[lag1]+_b[lag2]+_b[lag3]+_b[lag4]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
ereturn display
end


/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/************************************ ESTIMATION PROCEDURES AND STORING RESULTS ************************************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/

use ${project}/DDCGdata_final.dta, clear

/* Within estimates */
xtreg y l(1/4).y dem yy*   , fe r cluster(wbcode2)
estimates store e1
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e1_add


/* Abond estimates */
xtabond2 y l(1/4).y dem yy*   , gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*   , p) noleveleq robust nodiffsargan
estimates store e2
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e2_add

/* Abond with regular instruments for democracy*/
xtabond2 y l(1/4).y dem yy*   , gmmstyle(y, laglimits(2 .)) ivstyle(l.dem yy*   , p) noleveleq robust nodiffsargan
estimates store e3
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e3_add

/* Abond with truncated instruments for GDP*/
xtabond2 y l(1/4).y dem yy*   , gmmstyle(y, laglimits(2 25)) ivstyle(l.dem yy*   , p) noleveleq robust nodiffsargan
estimates store e4
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e4_add

/* Abond with truncated instruments for GDP in forward orthogonal differences */
xtabond2 y l(1/4).y dem yy*   , gmmstyle(y, laglimits(2 5)) ivstyle(yy*   l.dem, p) noleveleq robust orth nodiffsargan
estimates store e5
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e5_add

/* Ahn and Schmidt moments */
use ${project}/DDCGdata_final.dta, clear
local iter=10 /*Warning, this choice may affect the AS estimates since they may cicle and not converge*/
gen temp=.
replace temp=year if y!=.
bysort wbcode: egen maxy=max(temp)
gen temp0y=.
gen temp1y=.
gen temp2y=.
gen temp3y=.
gen temp4y=.
gen temp0d=.
replace temp0y=y if year==maxy
replace temp1y=y if year==maxy-1
replace temp2y=y if year==maxy-2
replace temp3y=y if year==maxy-3
replace temp4y=y if year==maxy-4
replace temp0d=dem if year==maxy
bysort wbcode2: egen y0=total(temp0y)
bysort wbcode2: egen y1=total(temp1y)
bysort wbcode2: egen y2=total(temp2y)
bysort wbcode2: egen y3=total(temp3y)
bysort wbcode2: egen y4=total(temp4y)
bysort wbcode2: egen d0=total(temp0d)

sort wbcode2 year
/* Iterative procedure to add Ahn and Schmidt moments which improve the estimator near a unit root */
xtreg y l(1/4).y dem yy*   , fe r cluster(wbcode2)
gen instrument=y0-_b[L.y]*y1-_b[L2.y]*y2-_b[L3.y]*y3-_b[L4.y]*y4-_b[dem]*d0
forvalues j=1(1)`iter'{
xtabond2 y l(1/4).y dem yy*   , gmmstyle(instrument, laglimits(1 1)) gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*   , p) noleveleq robust nodiffsargan
replace instrument=y0-_b[L.y]*y1-_b[L2.y]*y2-_b[L3.y]*y3-_b[L4.y]*y4-_b[dem]*d0
}
estimates store e6
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e6_add


/* Iterative procedure to add Ahn and Schmidt moments which improve the estimator near the unit circle */
xtreg y l(1/4).y dem yy*   , fe r cluster(wbcode2)
replace instrument=y0-_b[L.y]*y1-_b[L2.y]*y2-_b[L3.y]*y3-_b[L4.y]*y4-_b[dem]*d0
forvalues j=1(1)`iter'{
xtabond2 y l(1/4).y dem yy*   , gmmstyle(instrument, laglimits(1 1)) gmmstyle(y, laglimits(2 .)) ivstyle(yy*   l.dem, p) noleveleq robust nodiffsargan
replace instrument=y0-_b[L.y]*y1-_b[L2.y]*y2-_b[L3.y]*y3-_b[L4.y]*y4-_b[dem]*d0
}
estimates store e7
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e7_add

/* Iterative procedure to add Ahn and Schmidt moments which improve the estimator near the unit circle */
xtreg y l(1/4).y dem yy*   , fe r cluster(wbcode2)
replace instrument=y0-_b[L.y]*y1-_b[L2.y]*y2-_b[L3.y]*y3-_b[L4.y]*y4-_b[dem]*d0
forvalues j=1(1)`iter'{
xtabond2 y l(1/4).y dem yy*   , gmmstyle(instrument, laglimits(1 1)) gmmstyle(y, laglimits(2 25)) ivstyle(yy*   l.dem, p) noleveleq robust nodiffsargan
replace instrument=y0-_b[L.y]*y1-_b[L2.y]*y2-_b[L3.y]*y3-_b[L4.y]*y4-_b[dem]*d0
}
estimates store e8
nlcom (shortrun: _b[dem])  (lag1: _b[L1.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e8_add

/* Export table */
estout e1 e2 e3 e4 e5 e6 e7 e8  using  ${project}/appendix/TableGMM.tex, style(tex) ///
varlabels(L.y "log GDP first lag" L2.y "log GDP second lag" L3.y "log GDP third lag" L4.y "log GDP fourth lag"  dem "Democracy") ///
cells(b(star fmt(%9.3f)) se(par)) stats( ar2p j N N_g, fmt(%7.2f %7.0f %7.0f %7.0f) ///
labels("\input{appendix/TableGMM_Add}  AR2 test p-value" Moments Observations "Countries in sample"))   collabels(none) ///
keep(L.y L2.y L3.y L4.y dem) order(dem L.y L2.y L3.y L4.y) stardrop(dem L.y L2.y L3.y L4.y) nolabel replace mlabels(none)

estout e1_add e2_add e3_add e4_add e5_add e6_add e7_add e8_add using ${project}/appendix/TableGMM_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") ///
cells(b(star fmt(%9.3f)) se(par)) ///
keep(longrun effect$limit persistence) order(longrun effect$limit persistence) ///
stardrop(longrun effect$limit persistence)  nolabel replace mlabels(none) collabels(none) ///


