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

/* Define lags */
gen L1y=l.y
gen L2y=l2.y
gen L3y=l3.y
gen L4y=l4.y

/* Calculate cook's distance */
reg y dem L1y L2y L3y L4y  yy* i.wbcode2
predict cd, cooksd
predict er, rstandard
local num=e(N)

/* Calculate robust regression weights */
rreg y dem L1y L2y L3y L4y  yy*   i.wbcode2, gen(wr)


/* Estimates */
xtreg y L1y L2y L3y L4y dem yy*  , fe r cluster(wbcode2)
gen samp=e(sample)
estimates store e1
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e1_add

xtreg y dem L1y L2y L3y L4y  yy*  if abs(er)<1.96, fe cluster(wbcode2) r
estimates store e2
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e2_add

xtreg y dem L1y L2y L3y L4y  yy*  if cd<4/`num', fe cluster(wbcode2) r
estimates store e3
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e3_add

reg y dem L1y L2y L3y L4y  yy*   i.wbcode2 [aw=wr], cluster(wbcode2) r
estimates store e4
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e4_add

xi: mregress y dem L1y L2y L3y L4y  i.year  i.wbcode2 
estimates store e5
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e5_add

/* Exports tables */
estout e1 e2 e3 e4 e5 using ${project}/appendix/TableOut.tex, style(tex) ///
varlabels(L1y "log GDP first lag" L2y "log GDP second lag" L3y "log GDP third lag" L4y "log GDP fourth lag"  dem "Democracy") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N , fmt(%7.0f) labels( "\input{appendix/TableOut_Add} Observations"))  ///
collabels(none)  keep(dem L1y L2y L3y L4y ) order(dem L1y L2y L3y L4y ) stardrop(dem L1y L2y L3y L4y )   ///
nolabel replace mlabels(none)

estout e1_add e2_add e3_add e4_add e5_add using ${project}/appendix/TableOut_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") ///
cells(b(star fmt(%9.3f)) se(par)) ///
keep(longrun effect$limit persistence) ///
order(longrun effect$limit persistence) ///
stardrop(longrun effect$limit persistence) /// 
nolabel replace mlabels(none) collabels(none) /// 

