clear all
set matsize 5000
set more off
cap ssc install parmest
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */
local reps=50

use ${project}/DDCGdata_final.dta

gen trever=.
replace trever=1 if dem==0&l.dem==1
replace trever=0 if dem==1&l.dem==1

gen lag1y=l1.y
gen lag2y=l2.y
gen lag3y=l3.y
gen lag4y=l4.y

keep if lag1y!=.
keep if lag2y!=.
keep if lag3y!=.
keep if lag4y!=.

tab year, gen(dyear)

gen ydep0=L15.y-L.y
gen ydep1=L14.y-L.y
gen ydep2=L13.y-L.y
gen ydep3=L12.y-L.y
gen ydep4=L11.y-L.y
gen ydep5=L10.y-L.y
gen ydep6=L9.y-L.y
gen ydep7=L8.y-L.y
gen ydep8=L7.y-L.y
gen ydep9=L6.y-L.y
gen ydep10=L5.y-L.y
gen ydep11=L4.y-L.y
gen ydep12=L3.y-L.y
gen ydep13=L2.y-L.y
gen ydep14=0
gen ydep15=y-L.y
gen ydep16=F1.y-L.y
gen ydep17=F2.y-L.y
gen ydep18=F3.y-L.y
gen ydep19=F4.y-L.y
gen ydep20=F5.y-L.y
gen ydep21=F6.y-L.y
gen ydep22=F7.y-L.y
gen ydep23=F8.y-L.y
gen ydep24=F9.y-L.y
gen ydep25=F10.y-L.y
gen ydep26=F11.y-L.y
gen ydep27=F12.y-L.y
gen ydep28=F13.y-L.y
gen ydep29=F14.y-L.y
gen ydep30=F15.y-L.y
gen ydep31=F16.y-L.y
gen ydep32=F17.y-L.y
gen ydep33=F18.y-L.y
gen ydep34=F19.y-L.y
gen ydep35=F20.y-L.y
gen ydep36=F21.y-L.y
gen ydep37=F22.y-L.y
gen ydep38=F23.y-L.y
gen ydep39=F24.y-L.y
gen ydep40=F25.y-L.y
gen ydep41=F26.y-L.y
gen ydep42=F27.y-L.y
gen ydep43=F28.y-L.y
gen ydep44=F29.y-L.y
gen ydep45=F30.y-L.y

keep if trever!=.
set seed 12345


/**********************************************Regression adjustment results*****************************************************************/
program tefpas_ra, eclass

quietly: tab year if ydep0!=. & trever!=., gen(dumyears)

cap: teffects ra (ydep0 lag1y lag2y lag3y lag4y dumyears*, noconstant) (trever), atet iterate(7)
if _rc==0{
matrix b=_b[ATET: r1vs0.trever]
}
else{
matrix b=(.)
}

quietly: drop dumyears*

forvalues s=1(1)10{

quietly: tab year if ydep`s'!=. & trever!=., gen(dumyears)

cap: teffects ra (ydep`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (trever), atet iterate(7)
if _rc==0{
matrix b=(b, _b[ATET: r1vs0.trever])
}
else{
matrix b=(b, .)
}

quietly: drop dumyears*
}

matrix b=(b, 0)
matrix b=(b, 0)
matrix b=(b, 0)
matrix b=(b, 0)

forvalues s=15(1)45{

quietly: tab year if ydep`s'!=. & trever!=., gen(dumyears)

cap: teffects ra (ydep`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (trever), atet iterate(7)
if _rc==0{
matrix b=(b, _b[ATET: r1vs0.trever])
}
else{
matrix b=(b, .)
}

quietly: drop dumyears*
}

ereturn post b
end

tefpas_ra
xtset, clear
bootstrap _b, reps(`reps') cluster(wbcode2): tefpas_ra 
parmest, format(estimate min95 max95) saving(${project}/appendix/impulse_rever_ra, replace)



/**************************************Inverse Propensity Score Weighting Results*****************************************************************/
program tefpas_ipw, eclass

quietly: gen temp=trever if ydep0!=. &trever!=.
quietly: bysort year: egen mtemp=max(temp)
quietly: tab year if mtemp==1, gen(dumyears)

cap: teffects ipw (ydep0) (trever lag1y lag2y lag3y lag4y dumyears*, noconstant probit), atet iterate(7)
if _rc==0{
matrix b=_b[ATET: r1vs0.trever]
}
else{
matrix b=(.)
}

quietly: drop dumyears* temp mtemp

forvalues s=1(1)13{

quietly: gen temp=trever if ydep`s'!=. &trever!=.
quietly: bysort year: egen mtemp=max(temp)
quietly: tab year if mtemp==1, gen(dumyears)

cap: teffects ipw (ydep`s') (trever lag1y lag2y lag3y lag4y dumyears*, noconstant probit), atet iterate(7)
if _rc==0{
matrix b=(b, _b[ATET: r1vs0.trever])
}
else{
matrix b=(b, .)
}

quietly: drop dumyears* temp mtemp

}

matrix b=(b, 0)

forvalues s=15(1)45{

quietly: gen temp=trever if ydep`s'!=. &trever!=.
quietly: bysort year: egen mtemp=max(temp)
quietly: tab year if mtemp==1, gen(dumyears)

cap: teffects ipw (ydep`s') (trever lag1y lag2y lag3y lag4y dumyears*, noconstant probit), atet iterate(7)
if _rc==0{
matrix b=(b, _b[ATET: r1vs0.trever])
}
else{
matrix b=(b, .)
}

quietly: drop dumyears* temp mtemp

}

ereturn post b
end


tefpas_ipw
xtset, clear
bootstrap _b, reps(`reps') cluster(wbcode2): tefpas_ipw
parmest, format(estimate min95 max95) saving(${project}/appendix/impulse_rever_ipw, replace)



/*************************************************Regression adjustment and IPW results*********************************************************************/
program tefpas_ipwra, eclass

quietly: gen temp=trever if ydep0!=. &trever!=.
quietly: bysort year: egen mtemp=max(temp)
quietly: tab year if mtemp==1, gen(dumyears)

cap: teffects ipwra (ydep0 lag1y lag2y lag3y lag4y dumyears*, noconstant) (trever lag1y lag2y lag3y lag4y dumyears*, probit), atet iterate(7)
if _rc==0{
matrix b=_b[ATET: r1vs0.trever]
}
else{
matrix b=(.)
}

quietly: drop dumyears* temp mtemp

forvalues s=1(1)10{

quietly: gen temp=trever if ydep`s'!=. &trever!=.
quietly: bysort year: egen mtemp=max(temp)
quietly: tab year if mtemp==1, gen(dumyears)

cap: teffects ipwra (ydep`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (trever lag1y lag2y lag3y lag4y dumyears*, probit), atet iterate(7)
if _rc==0{
matrix b=(b, _b[ATET: r1vs0.trever])
}
else{
matrix b=(b, .)
}

quietly: drop dumyears* temp mtemp
}

matrix b=(b, 0)
matrix b=(b, 0)
matrix b=(b, 0)
matrix b=(b, 0)

forvalues s=15(1)45{

quietly: gen temp=trever if ydep`s'!=. &trever!=.
quietly: bysort year: egen mtemp=max(temp)
quietly: tab year if mtemp==1, gen(dumyears)

cap: teffects ipwra (ydep`s' lag1y lag2y lag3y lag4y dumyears*, noconstant) (trever lag1y lag2y lag3y lag4y dumyears*, probit), atet iterate(7)
if _rc==0{
matrix b=(b, _b[ATET: r1vs0.trever])
}
else{
matrix b=(b, .)
}

quietly: drop dumyears* temp mtemp
}

ereturn post b
end


xtset, clear
bootstrap _b, reps(`reps') cluster(wbcode2): tefpas_ipwra
parmest, format(estimate min95 max95) saving(${project}/appendix/impulse_rever_ipwra, replace)


/*******************************************************Figure: local projections RA**********************************************************/
use ${project}/appendix/impulse_rever_ra, clear
split parm, p("c")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Change in GDP per capita log points") xtitle(Years around reversal) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)30) ylabel(-40(10)20)
graph export ${project}/appendix/ImpRevRA.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/

/*Figure: local projections IPW*/
use ${project}/appendix/impulse_rever_ipw, clear
split parm, p("c")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Change in GDP per capita log points") xtitle(Years around reversal) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)30) ylabel(-40(10)20)
graph export ${project}/appendix/ImpRevIPW.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/


/*Figure: local projections IPWRA*/
use ${project}/appendix/impulse_rever_ipwra, clear
split parm, p("c")
keep if parm2!=""
destring parm2, force replace
gen time=parm2-16
tsset time

twoway (tsline estimate, lcolor(black) lpattern(solid)) ///
       (tsline min95, lcolor(gray) lpattern(dash)) ///
       (tsline max95, lcolor(gray) lpattern(dash)) ///
       , ytitle("Change in GDP per capita log points") xtitle(Years around reversal) legend(off) ///
	   graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ///
	   xlabel(-15(5)30) ylabel(-40(10)20)
graph export ${project}/appendix/ImpRevIPWRA.eps, as(eps) preview(off) replace /*Plots and saves figure 1*/


