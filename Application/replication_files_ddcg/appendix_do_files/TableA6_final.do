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


/* Loop over different measures of democracy */
foreach var in dem demPS demFH demPOL demCGV demBMR{

/*Construct instrument (waves) for for alternative democracy measure*/
gen demest=`var'
gen tempy=.
replace tempy=year if demest!=.
bysort wbcode: egen fy=min(tempy)
gen temp=.
replace temp=1 if demest==1&year==fy
replace temp=1 if demest==1&year==fy+1
replace temp=1 if demest==1&year==fy+2
replace temp=1 if demest==1&year==fy+3
replace temp=1 if demest==1&year==fy+4
bysort wbcode: egen tagtemp=total(temp)
gen initdemFYdemest=0
replace initdemFYdemest=1 if tagtemp==5
drop temp tagtemp fy tempy

gen demextdemest="_dem"
replace demextdemest="_nd" if initdemFYdemest==0
sort wbcode2 year
gen region2demest=region+demextdemest

bysort region2demest year: egen total=total(demest)
gen jtotal=total-demest
replace jtotal=total if demest==.
bysort region2demest year: egen count=count(demest)
gen jcount=count-1
replace jcount=count if demest==.
gen demestreg=jtotal/jcount
drop total count jtotal jcount

sort wbcode2 year 

xtreg   y  demest yy*  , fe r cluster(wbcode2)
estimates store epan_`var'

xtreg   y  l(1/4).y demest yy*  , fe r cluster(wbcode2)
estimates store edyn_`var'
nlcom (shortrun: _b[demest])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store edyn_`var'_add

xtivreg2 y l(1/4).y yy* (demest=l(1/4).demestreg), fe cluster(wbcode2) r partial(yy*)
estimates store eiv_`var'
nlcom (shortrun: _b[demest])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
vareffects
estimates store eiv_`var'_add

drop demest  demestreg initdemFYdemest region2demest demextdemest

}

/* Export tables */
#delimit ;
estout edyn_dem edyn_demPS edyn_demFH edyn_demPOL edyn_demCGV edyn_demBMR using ${project}/appendix/TableOthersTop.tex, style(tex) 
varlabels(demest "Democracy") 
cells(b(star fmt(%9.3f)) se(par)) 
stats(N N_g, fmt(%7.0f  %7.0f) labels( "\input{appendix/TableOthersTop_Add} Observations" "Countries in sample" ))  
keep(demest) 
order(demest)
stardrop(demest)  
nolabel replace mlabels(none) collabels(none); 
#delimit cr

#delimit ;
estout edyn_dem_add edyn_demPS_add edyn_demFH_add edyn_demPOL_add edyn_demCGV_add edyn_demBMR_add using ${project}/appendix/TableOthersTop_Add.tex, style(tex) 
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence) 
order(longrun effect$limit persistence)
stardrop(longrun effect$limit persistence)  
nolabel replace mlabels(none) collabels(none); 
#delimit cr

#delimit ;
estout eiv_dem eiv_demPS eiv_demFH eiv_demPOL eiv_demCGV eiv_demBMR using ${project}/appendix/TableOthersMid.tex, style(tex) 
varlabels(demest "Democracy") 
cells(b(star fmt(%9.3f)) se(par)) 
stats(N N_g, fmt(%7.0f  %7.0f) labels( "\input{appendix/TableOthersMid_Add} Observations" "Countries in sample" ))  
keep(demest) 
order(demest)
stardrop(demest)  
nolabel replace mlabels(none) collabels(none); 
#delimit cr

#delimit ;
estout eiv_dem_add eiv_demPS_add eiv_demFH_add eiv_demPOL_add eiv_demCGV_add eiv_demBMR_add using ${project}/appendix/TableOthersMid_Add.tex, style(tex) 
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence) 
order(longrun effect$limit persistence)
stardrop(longrun effect$limit persistence)  
nolabel replace mlabels(none) collabels(none); 
#delimit cr


#delimit ;
estout epan_dem epan_demPS epan_demFH epan_demPOL epan_demCGV epan_demBMR using ${project}/appendix/TableOthersBot.tex, style(tex) 
varlabels(demest "Democracy") 
cells(b(star fmt(%9.3f)) se(par)) 
stats(N N_g, fmt(%7.0f  %7.0f) labels( " Observations" "Countries in sample" ))  
keep(demest) 
order(demest)
stardrop(demest)  
nolabel replace mlabels(none) collabels(none); 
#delimit cr
