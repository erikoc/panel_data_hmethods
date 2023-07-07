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

use ${project}/DDCGdata_final.dta

gen L1y=l.y
gen L2y=l2.y
gen L3y=l3.y
gen L4y=l4.y
gen ldemreg1=l.demreg
gen ldemreg2=l2.demreg
gen ldemreg3=l3.demreg
gen ldemreg4=l4.demreg
tab wbcode, gen(www)

/*Calculate prediction errors*/
ivreg2 y l(1/4).y (dem=l(1/4).demreg) yy*   www*, cluster(wbcode2) r
gen sample=e(sample)
predict er2, stdp

/*Calculate cook's distance*/
reg dem l(1/4).y L(1/4).demreg yy*   i.wbcode2 if sample==1
predict predicted_dem, xb

reg y l(1/4).y yy*   predicted_dem i.wbcode2 if sample==1
predict cd, cooksd
predict er, rstandard
local num=e(N)

/*Calculate robust regression weights*/
rreg y l(1/4).y yy*   predicted_dem i.wbcode2 if sample==1, gen(wr)

/*****************************************************/
/**Estimates removing outliers from the second stage**/
/*****************************************************/
xtivreg2 y L1y L2y L3y L4y (dem=l(1/4).demreg) yy*  , fe r cluster(wbcode2)
estimates store e1
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e1_add


xtivreg2 y L1y L2y L3y L4y (dem=L(1/4).demreg) yy*   if abs(er)<1.96&sample==1, fe cluster(wbcode2) r
estimates store e2
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e2_add

xtivreg2 y L1y L2y L3y L4y (dem=L(1/4).demreg) yy*   if cd<4/`num'&sample==1, fe cluster(wbcode2) r
estimates store e3
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e3_add


xtivreg2 y L1y L2y L3y L4y (dem=L(1/4).demreg) yy*    if sample==1 [aw=wr],fe cluster(wbcode2) r
estimates store e4
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e4_add


/********************************************************************/
/**Estimates removing outliers from both the first and second stage**/
/**Includes programs to adjust the variance of these estimators *****/
/********************************************************************/

rename dem z
gen cons=1
keep if sample==1

/* Model 1 */
capture: drop dem 
capture: drop wnew*
capture: drop wnewa*
capture: drop zres
capture: drop sigma2z
capture: drop zz
capture: drop yhat
capture: drop yres
capture: drop sigma2y
capture: drop dd

reg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
predict sterr, rstand
tab wbcode if e(sample)==1&abs(sterr)<1.96, gen(wnew)

reg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 if sterr<1.96, noconstant cluster(wbcode2) r
predict dem, xb
gen zres=z-dem
sum zres
gen sigma2z=r(Var)
matrix V1=e(V)

reg y dem L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
predict sterr2, rstand 
tab wbcode if e(sample)==1&abs(sterr2)<1.96, gen(wnewa)

reg y dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173 if abs(sterr2)<1.96, noconstant cluster(wbcode2) r
matrix V2=e(V)
scalar zz=_b[dem]
predict yhat, xb
gen yres=y-yhat
sum yres
gen sigma2y=r(Var)

/*Two-step adjusted variance matrix*/
matrix accum R =  ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173  [iw=(y-yhat)*(z-dem)/(sigma2z*sigma2y)], nocons
matrix R=R[229..453,1..228]
matrix accum Cs1 = ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173 [iw=-zz/sigma2y], nocons
matrix Cs1 = Cs1[229..453,1..228] /* Get only the desired partition */
gen dd = (y-yhat)/sigma2y
matrix vecaccum Cs2 = dd ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173, nocons
matrix Cs2 = Cs2' , J(228,224,0) /* Plug into the relevant column */
matrix Cs = -(Cs1 + Cs2')
matrix Ms = V2 + (V2 * (Cs*V1*Cs' - R*V1*Cs' - Cs*V1*R') * V2)
matrix b = e(b)
capture program drop poster
program poster, eclass
ereturn post b Ms /* For sandwich results: est post b Ms */
ereturn local vcetype "Huber M"
ereturn local cmd "Robust IV"
end
/*Display estimates*/
poster
ereturn display
estimates store e5
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e5_add



/* Model 2 */
capture: drop dem 
capture: drop wnew*
capture: drop wnewa*
capture: drop zres
capture: drop sigma2z
capture: drop zz
capture: drop yhat
capture: drop yres
capture: drop sigma2y
capture: drop dd

reg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
predict cdist, cooks 
local num=e(N)
tab wbcode if e(sample)==1&cdist<4/`num', gen(wnew)

reg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 if cdist<4/`num', noconstant cluster(wbcode2) r
predict dem, xb
gen zres=z-dem
sum zres
gen sigma2z=r(Var)
matrix V1=e(V)

reg y dem L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
predict cdist2, cooks 
local num=e(N)
tab wbcode if e(sample)==1&cdist2<4/`num', gen(wnewa)

reg y dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173 if cdist2<4/`num', noconstant cluster(wbcode2) r
matrix V2=e(V)
scalar zz=_b[dem]
predict yhat, xb
gen yres=y-yhat
sum yres
gen sigma2y=r(Var)

/*Two-step adjusted variance matrix*/
matrix accum R =  ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173  [iw=(y-yhat)*(z-dem)/(sigma2z*sigma2y)], nocons
matrix R=R[229..453,1..228]
matrix accum Cs1 = ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173 [iw=-zz/sigma2y], nocons
matrix Cs1 = Cs1[229..453,1..228] /* Get only the desired partition */
gen dd = (y-yhat)/sigma2y
matrix vecaccum Cs2 = dd ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173, nocons
matrix Cs2 = Cs2' , J(228,224,0) /* Plug into the relevant column */
matrix Cs = -(Cs1 + Cs2')
matrix Ms = V2 + (V2 * (Cs*V1*Cs' - R*V1*Cs' - Cs*V1*R') * V2)
matrix b = e(b)
capture program drop poster
program poster, eclass
ereturn post b Ms /* For sandwich results: est post b Ms */
ereturn local vcetype "Huber M"
ereturn local cmd "Robust IV"
end
/*Display estimates*/
poster
ereturn display
estimates store e6
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e6_add




/* Model 3 */
capture: drop dem 
capture: drop wnew*
capture: drop wnewa*
capture: drop zres
capture: drop sigma2z
capture: drop zz
capture: drop yhat
capture: drop yres
capture: drop sigma2y
capture: drop dd

reg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
tab wbcode if e(sample), gen(wnew)

rreg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y yy6-yy51 wnew1-wnew173
predict dem, xb
gen zres=z-dem
sum zres
gen sigma2z=r(Var)
matrix V1=e(V)

reg y dem L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
tab wbcode if e(sample), gen(wnewa)

rreg y dem L1y L2y L3y L4y yy6-yy51 wnewa1-wnewa173
matrix V2=e(V)
scalar zz=_b[dem]
predict yhat, xb
gen yres=y-yhat
sum yres
gen sigma2y=r(Var)

/*Two-step adjusted variance matrix*/
matrix accum R =  ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173  [iw=(y-yhat)*(z-dem)/(sigma2z*sigma2y)], nocons
matrix R=R[229..453,1..228]
matrix accum Cs1 = ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173 [iw=-zz/sigma2y], nocons
matrix Cs1 = Cs1[229..453,1..228] /* Get only the desired partition */
gen dd = (y-yhat)/sigma2y
matrix vecaccum Cs2 = dd ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173, nocons
matrix Cs2 = Cs2' , J(228,224,0) /* Plug into the relevant column */
matrix Cs = -(Cs1 + Cs2')
matrix Ms = V2 + (V2 * (Cs*V1*Cs' - R*V1*Cs' - Cs*V1*R') * V2)
matrix b = e(b)
capture program drop poster
program poster, eclass
ereturn post b Ms /* For sandwich results: est post b Ms */
ereturn local vcetype "Huber M"
ereturn local cmd "Robust IV"
end
/*Display estimates*/
poster
ereturn display
estimates store e7
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e7_add


/* Model 4 */
capture: drop dem 
capture: drop wnew*
capture: drop wnewa*
capture: drop zres
capture: drop sigma2z
capture: drop zz
capture: drop yhat
capture: drop yres
capture: drop sigma2y
capture: drop dd

reg z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
tab wbcode if e(sample), gen(wnew)

xi: mregress z ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173, noconstant
predict dem, xb
gen zres=z-dem
sum zres
gen sigma2z=r(Var)
matrix V1=e(V)

reg y dem L1y L2y L3y L4y cons yy6-yy51 www*, noconstant
tab wbcode if e(sample), gen(wnewa)

xi: mregress y dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173, noconstant
matrix V2=e(V)
scalar zz=_b[dem]
predict yhat, xb
gen yres=y-yhat
sum yres
gen sigma2y=r(Var)

/*Two-step adjusted variance matrix*/
matrix accum R =  ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173  [iw=(y-yhat)*(z-dem)/(sigma2z*sigma2y)], nocons
matrix R=R[229..453,1..228]
matrix accum Cs1 = ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173 dem L1y L2y L3y L4y cons yy6-yy51 wnewa1-wnewa173 [iw=-zz/sigma2y], nocons
matrix Cs1 = Cs1[229..453,1..228] /* Get only the desired partition */
gen dd = (y-yhat)/sigma2y
matrix vecaccum Cs2 = dd ldemreg1 ldemreg2 ldemreg3 ldemreg4 L1y L2y L3y L4y cons yy6-yy51 wnew1-wnew173, nocons
matrix Cs2 = Cs2' , J(228,224,0) /* Plug into the relevant column */
matrix Cs = -(Cs1 + Cs2')
matrix Ms = V2 + (V2 * (Cs*V1*Cs' - R*V1*Cs' - Cs*V1*R') * V2)
matrix b = e(b)
capture program drop poster
program poster, eclass
ereturn post b Ms /* For sandwich results: est post b Ms */
ereturn local vcetype "Huber M"
ereturn local cmd "Robust IV"
end
/*Display estimates*/
poster
ereturn display
estimates store e8
nlcom (shortrun: _b[dem])  (lag1: _b[L1y])  (lag2: _b[L2y])  (lag3: _b[L3y])  (lag4: _b[L4y]), post
vareffects
estimates store e8_add

/* Export tables */
estout e1 e2 e3 e4 e5 e6 e7 e8 using ${project}/appendix/TableOutliers2SLS.tex, style(tex) ///
varlabels(L1y "log GDP first lag" L2y "log GDP second lag" L3y "log GDP third lag" L4y "log GDP fourth lag"  dem "Democracy") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N , fmt(%7.0f) labels( "\input{appendix/TableOutliers2SLS_Add} Observations"))  collabels(none)  ///
keep(L1y L2y L3y L4y dem) order(dem L1y L2y L3y L4y) stardrop(dem L1y L2y L3y L4y)   nolabel replace mlabels(none)

estout e1_add e2_add e3_add e4_add e5_add e6_add e7_add e8_add  using ${project}/appendix/TableOutliers2SLS_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") /// 
cells(b(star fmt(%9.3f)) se(par)) keep(longrun effect$limit persistence) order(longrun effect$limit persistence) stardrop(longrun effect$limit persistence)  ///
nolabel replace mlabels(none) collabels(none)
