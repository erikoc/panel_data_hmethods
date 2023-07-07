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

/* Determines sample */
xtreg y l(1/4).y dem yy*, fe
gen sampOLS=e(sample)

/* Interactions with the share of population with primary schooling */
gen temp=lp_bl+ls_bl+lh_bl if year==1960
bysort wbcode: egen sec1960=mean(temp)
gen var=sec1960 if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter  yy*, fe r cluster(wbcode2)
estimates store e1A
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e1A_add
drop var  inter temp sec1960

gen temp=lp_bl+ls_bl+lh_bl if year==1970
bysort wbcode: egen sec1970=mean(temp)
gen var=sec1970 if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter  yy*, fe r cluster(wbcode2)
estimates store e2A
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e2A_add
drop var  inter temp sec1970

gen temp=lp_bl+ls_bl+lh_bl if year==1980
bysort wbcode: egen sec1980=mean(temp)
gen var=sec1980 if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter  yy*, fe r cluster(wbcode2)
estimates store e3A
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e3A_add
drop var  inter temp sec1980

gen temp=lp_bl+ls_bl+lh_bl
gen quinq=floor(year/5) 
bysort wbcode quinq: egen secondary=mean(temp)
sort wbcode2 year
gen var=l.secondary if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter var yy*, fe r cluster(wbcode2)
estimates store e4A
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e4A_add
drop var  inter temp secondary quinq 

/* Interactions with the share of population with tertiary education */
gen temp=lh_bl if year==1960
bysort wbcode: egen sec1960=mean(temp)
gen var=sec1960 if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter  yy*, fe r cluster(wbcode2)
estimates store e1B
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e1B_add
drop var  inter temp sec1960

gen temp=lh_bl if year==1970
bysort wbcode: egen sec1970=mean(temp)
gen var=sec1970 if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter  yy*, fe r cluster(wbcode2)
estimates store e2B
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e2B_add
drop var  inter temp sec1970

gen temp=lh_bl if year==1980
bysort wbcode: egen sec1980=mean(temp)
gen var=sec1980 if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter  yy*, fe r cluster(wbcode2)
estimates store e3B
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e3B_add
drop var  inter temp sec1980

gen temp=lh_bl
gen quinq=floor(year/5) 
bysort wbcode quinq: egen secondary=mean(temp)
sort wbcode2 year
gen var=l.secondary if sampOLS==1
sum var, d
local varev=r(p25)
gen inter=dem*(var-`varev')
sort wbcode2 year

xtreg y l(1/4).y dem inter var yy*, fe r cluster(wbcode2)
estimates store e4B
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e4B_add
drop var  inter temp secondary quinq 

/* Export tables */
estout e1A e2A e3A e4A e1B e2B e3B e4B   using ${project}/appendix/TableInterEduc.tex, style(tex) varlabels(dem "Democracy" inter "Interaction") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N N_g, fmt(%7.0f %7.0f) labels("\input{appendix/TableInterEduc_Add} Observations" "Countries in sample" )) ///
 collabels(none)  keep(dem inter) order(dem inter)  nolabel replace mlabels(none)  stardrop(dem inter)

estout e1A_add e2A_add e3A_add e4A_add e1B_add e2B_add e3B_add e4B_add  using ${project}/appendix/TableInterEduc_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") ///
cells(b(star fmt(%9.3f)) se(par)) keep(longrun effect$limit persistence) order(longrun effect$limit persistence) stardrop(longrun effect$limit persistence) ///
nolabel replace mlabels(none) collabels(none)


