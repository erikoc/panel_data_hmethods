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

/*Number of cells initial regime/region*/
local n60=13 /*initial regime in 1960*/
local nDA=12 /*Always democratic*/
local nREG=35 /*detailed regimes in 1960*/
local nFY=13 /*initial regime after independence*/

/*Bartik style instruments*/
bysort region year: egen total=total(dem)
gen jtotal=total-dem
replace jtotal=total if dem==.
bysort region year: egen count=count(dem)
gen jcount=count-1
replace jcount=count if dem==.
gen dembartikTemp=jtotal/jcount

forvalues j=1(1)`nFY'{
gen dembartikFY`j'=0
replace dembartikFY`j'=dembartikTemp if dFY_`j'==1
}


forvalues j=1(1)`n60'{
gen dembartik60`j'=0
replace dembartik60`j'=dembartikTemp if d60_`j'==1
}

forvalues j=1(1)`nDA'{
gen dembartikDA`j'=0
replace dembartikDA`j'=dembartikTemp if dDA_`j'==1
}

forvalues j=1(1)`nREG'{
gen dembartikREG`j'=0
replace dembartikREG`j'=dembartikTemp if dREG_`j'==1
}

drop total count jtotal jcount 


/*Base exercise*/
xtivreg2 y l(1/4).y (dem=l(1/4).demreg) yy*, fe cluster(wbcode2) r partial(yy*)
gen sampleFY=e(sample)
estimates store e1
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e1_add


/*Initial regime using 1960-1964 info*/
xtivreg2 y l(1/4).y (dem=l(1/4).demreg60) yy*, fe cluster(wbcode2) r  partial(yy*)
gen sample60=e(sample)
estimates store e2
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e2_add



/*Initial regime using democratic always*/
xtivreg2 y l(1/4).y (dem=l(1/4).demregDA) yy*, fe cluster(wbcode2) r  partial(yy*)
gen sampleDA=e(sample)
estimates store e3
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e3_add


/*Regimes instrument*/
xtivreg2 y l(1/4).y (dem=l(1/4).demregREGIME) yy* , fe cluster(wbcode2) r  partial(yy*)
gen sampleREGIME=e(sample)
estimates store e4
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e4_add


/* Bartik-style instruments */
xtivreg2 y l(1/4).y (dem=l(1/4).(dembartikFY1 dembartikFY2 dembartikFY3 dembartikFY4 dembartikFY5 dembartikFY6 dembartikFY7 dembartikFY8 dembartikFY9 dembartikFY10 dembartikFY11 dembartikFY12 dembartikFY13))  yy* if sampleFY==1, fe cluster(wbcode2) r partial(yy*)
estimates store e1bartik
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e1bartik_add


xtivreg2 y l(1/4).y (dem=l(1/4).(dembartik601 dembartik602 dembartik603 dembartik604 dembartik605 dembartik606 dembartik607 dembartik608 dembartik609 dembartik6010 dembartik6011 dembartik6012 dembartik6013))  yy* if sample60==1, fe cluster(wbcode2) r partial(yy*) 
estimates store e2bartik
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e2bartik_add


xtivreg2 y l(1/4).y (dem=l(1/4).(dembartikDA1 dembartikDA2 dembartikDA3 dembartikDA4 dembartikDA5 dembartikDA6 dembartikDA7 dembartikDA8 dembartikDA9 dembartikDA10 dembartikDA11 dembartikDA12))  yy* if sampleDA==1, fe cluster(wbcode2) r partial(yy*) 
estimates store e3bartik
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e3bartik_add


xtivreg2 y l(1/4).y (dem=l(1/4).(dembartikREG1	dembartikREG2	dembartikREG3	dembartikREG4	dembartikREG5	dembartikREG6	dembartikREG7	dembartikREG8	dembartikREG9	dembartikREG10	dembartikREG11	dembartikREG12	dembartikREG13	dembartikREG14	dembartikREG15	dembartikREG16	dembartikREG17	dembartikREG18	dembartikREG19	dembartikREG20	dembartikREG21	dembartikREG22	dembartikREG23	dembartikREG24	dembartikREG25	dembartikREG26	dembartikREG27	dembartikREG28	dembartikREG29	dembartikREG30	dembartikREG31	dembartikREG32	dembartikREG33 dembartikREG34 dembartikREG35)) yy* if sampleREGIME==1, fe cluster(wbcode2) r partial(yy*) 
estimates store e4bartik
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store e4bartik_add

/* Export tables */
estout e1 e2 e3 e4 e1bartik e2bartik e3bartik e4bartik using ${project}/appendix/TableIVs.tex, style(tex) varlabels(dem " Democracy") ///
cells(b(star fmt(%9.3f)) se(par)) stats( cdf  N N_g , fmt(%7.1f %7.0f %7.0f) labels("\input{appendix/TableIVs_Add} Exc. instruments F-stat." Observations "Countries in sample" ))  ///
 collabels(none)  keep(dem) order(dem) stardrop(dem) nolabel replace mlabels(none)

estout e1_add e2_add e3_add e4_add e1bartik_add e2bartik_add e3bartik_add e4bartik_add using ${project}/appendix/TableIVs_Add.tex, style(tex) ///
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") ///
cells(b(star fmt(%9.3f)) se(par)) keep(longrun effect$limit persistence) order(longrun effect$limit persistence) stardrop(longrun effect$limit persistence)  ///
nolabel replace mlabels(none) collabels(none)

