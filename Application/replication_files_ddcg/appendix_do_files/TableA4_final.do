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

if 1-_b[lag1]-_b[lag2]-_b[lag3]-_b[lag4]==0 {
quietly: nlcom (effect$limit: _b[effect$limit]) ///
      (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
ereturn display
}
else {
quietly: nlcom (effect$limit: _b[effect$limit]) ///
      (longrun: _b[shortrun]/(1-_b[lag1]-_b[lag2]-_b[lag3]-_b[lag4])) ///
      (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
ereturn display
}


end

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/************************************ ESTIMATION PROCEDURES AND STORING RESULTS ************************************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
use ${project}/DDCGdata_final.dta

/*Imposes persistence level of 0.95 */
sort wbcode2 year
local rho=.95
gen qdy=y-`rho'*l.y

xtreg qdy dem l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e1OLS
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e1OLS_add

xtivreg2 qdy (dem=l(1/4).demreg) l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e1IV
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e1IV_add
drop qdy 

/*Imposes persistence level of 0.96 */
sort wbcode2 year
local rho=.96
gen qdy=y-`rho'*l.y

xtreg qdy dem l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e2OLS
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e2OLS_add

xtivreg2 qdy (dem=l(1/4).demreg) l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e2IV
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e2IV_add
drop qdy 

/*Imposes persistence level of 0.97 */
sort wbcode2 year
local rho=.97
gen qdy=y-`rho'*l.y

xtreg qdy dem l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e3OLS
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e3OLS_add

xtivreg2 qdy (dem=l(1/4).demreg) l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e3IV
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e3IV_add
drop qdy 

/*Imposes persistence level of 0.98 */
sort wbcode2 year
local rho=.98
gen qdy=y-`rho'*l.y

xtreg qdy dem l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e4OLS
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e4OLS_add

xtivreg2 qdy (dem=l(1/4).demreg) l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e4IV
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e4IV_add
drop qdy 

/*Imposes persistence level of 0.99 */
sort wbcode2 year
local rho=.99
gen qdy=y-`rho'*l.y

xtreg qdy dem l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e5OLS
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e5OLS_add

xtivreg2 qdy (dem=l(1/4).demreg) l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e5IV
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e5IV_add
drop qdy 

/*Imposes persistence level of 1 (unit root) */
sort wbcode2 year
local rho=1
gen qdy=y-`rho'*l.y

xtreg qdy dem l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e6OLS
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e6OLS_add

xtivreg2 qdy (dem=l(1/4).demreg) l(1/3).d.y yy*, fe r cluster(wbcode2)
estimates store e6IV
nlcom (shortrun: _b[dem])  (lag1: `rho'+_b[LD.y])  (lag2: _b[L2D.y]-_b[LD.y])  (lag3: _b[L3D.y]-_b[L2D.y])  (lag4: -_b[L3D.y]), post
vareffects
estimates store e6IV_add
drop qdy 

/* Export tables */
estout e1OLS e2OLS e3OLS e4OLS e5OLS e6OLS  using ${project}/appendix/TableStationaryTop.tex, style(tex) ///
varlabels(dem "Democracy") cells(b(star fmt(%9.3f)) se(par)) ///
stats(N N_g, fmt(%7.0f %7.0f) labels( "\input{appendix/TableStationaryTop_Add}  Observations" "Countries in sample" ))  ///
keep(dem) order(dem) stardrop(dem)  nolabel replace mlabels(none) collabels(none) 

estout e1OLS_add e2OLS_add e3OLS_add e4OLS_add e5OLS_add e6OLS_add using ${project}/appendix/TableStationaryTop_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy" effect$limit "Effect of democracy after $limit years" ) ///
cells(b(star fmt(%9.3f)) se(par)) keep(longrun effect$limit) order(longrun effect$limit) stardrop(longrun effect$limit)  ///
nolabel replace mlabels(none) collabels(none)

estout e1IV e2IV e3IV e4IV e5IV e6IV  using ${project}/appendix/TableStationaryBot.tex, style(tex) ///
varlabels(dem "Democracy") cells(b(star fmt(%9.3f)) se(par)) stats(N N_g widstat, fmt(%7.0f %7.0f  %7.2f) ///
labels( "\input{appendix/TableStationaryBot_Add}  Observations" "Countries in sample" "Exc. Instruments F-stat."))  ///
keep(dem) order(dem) stardrop(dem)  nolabel replace mlabels(none) collabels(none) 

estout e1IV_add e2IV_add e3IV_add e4IV_add e5IV_add e6IV_add using ${project}/appendix/TableStationaryBot_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy" effect$limit "Effect of democracy after $limit years" )  ///
cells(b(star fmt(%9.3f)) se(par)) keep(longrun effect$limit) order(longrun effect$limit) stardrop(longrun effect$limit)  ///
nolabel replace mlabels(none) collabels(none)
