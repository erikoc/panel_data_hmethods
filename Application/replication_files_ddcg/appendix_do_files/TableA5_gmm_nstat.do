clear all
set more off
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */


/******************/
/*Install programs*/
/******************/
cap ssc install xtabond2 
mata: mata set matafavor speed


/*Defines program to simulate data and perform estimation*/
program define mcsimul, eclass
syntax [, stationary(integer 1) dimT(integer 40) dimN(integer 160) period(integer 20) beta(real 1) ar1(real 0.9) ar2(real 0) ar3(real 0) ar4(real 0) dar1(real -0.9) dar2(real 1) dar3(real 0) dar4(real 0)  vardem(real 1) vargdp(real 1) gdpshock(real 1) demshock(real 1) gdpinit(real 1) deminit(real 1)]
clear

set obs `dimN'
gen id_country=_n
expand `dimT'
bys id_country: gen id_time=_n

gen double u_y=invnorm(uniform())*`vargdp' if id_time==1
bys id_country: replace u_y=u_y[_n-1] if id_time>1

gen double u_dem=invnorm(uniform())*`vardem' if id_time==1
bys id_country: replace u_dem=u_dem[_n-1] if id_time>1

/* stationary or non-stationary initial conditions */
if `stationary'==1 {
gen double y=(u_y+u_dem)/(1-`ar1'-`ar2'-`ar3'-`ar4'-`dar1'-`dar2'-`dar3'-`dar4') if id_time==1
gen double dem=u_dem+y*(`dar1'+`dar2'+`dar3'+`dar4') if id_time==1
}
else{
gen double y=invnorm(uniform())*`gdpinit' if id_time==1
gen double dem=invnorm(uniform())*`deminit' if id_time==1
}

bys id_country: replace dem=u_dem+y[_n-1]*`dar1'+invnorm(uniform())*`demshock' if id_time==2
bys id_country: replace y=u_y+dem*`beta'+y[_n-1]*`ar1'+invnorm(uniform())*`gdpshock' if id_time==2

bys id_country: replace dem=u_dem+y[_n-1]*`dar1'+y[_n-2]*`dar2'+invnorm(uniform())*`demshock' if id_time==3
bys id_country: replace y=u_y+dem*`beta'+y[_n-1]*`ar1'+y[_n-2]*`ar2'+invnorm(uniform())*`gdpshock' if id_time==3

bys id_country: replace dem=u_dem+y[_n-1]*`dar1'+y[_n-2]*`dar2'+y[_n-3]*`dar3'+invnorm(uniform())*`demshock' if id_time==4
bys id_country: replace y=u_y+dem*`beta'+y[_n-1]*`ar1'+y[_n-2]*`ar2'+y[_n-3]*`ar3'+invnorm(uniform())*`gdpshock' if id_time==4

forvalues j=5(1)`dimT'{
bys id_country: replace dem=u_dem+y[_n-1]*`dar1'+y[_n-2]*`dar2'+y[_n-3]*`dar3'+y[_n-4]*`dar4'+invnorm(uniform())*`demshock' if id_time==`j'
bys id_country: replace y=u_y+dem*`beta'+y[_n-1]*`ar1'+y[_n-2]*`ar2'+y[_n-3]*`ar3'+y[_n-4]*`ar4'+invnorm(uniform())*`gdpshock' if id_time==`j'
}

keep if id_time>`period'
replace id_time=id_time-`period'
xtset id_country id_time
quietly: tab id_time, gen(yy)

quietly: xtabond2 y l(1/4).y dem yy* ,  gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*, p) noleveleq robust nodiffsargan

nlcom (shortrun_t: _b[dem]/_se[dem]) (shortrun: _b[dem]) (longrun: _b[dem]/(1-_b[L.y]-_b[L2.y]-_b[L3.y]-_b[L4.y])) (persistence: _b[L.y]+_b[L2.y]+_b[L3.y]+_b[L4.y])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]), post
matrix b=e(b)
ereturn post b
end

/* Get estimated values and empirical properties of the data required for the simulation */
use ${project}/DDCGdata_final.dta

xtreg y i.year l(1/4).y dem , fe r cluster(wbcode2)
local true_gdpshock=e(sigma_e)
local true_vargdp=e(sigma_u)
local true_ar1=_b[L1.y]
local true_ar2=_b[L2.y]
local true_ar3=_b[L3.y]
local true_ar4=_b[L4.y]
local true_beta=_b[dem]
local true_beta_t=_b[dem]/_se[dem]
local true_longest=_b[dem]/(1-_b[L1.y]-_b[L2.y]-_b[L3.y]-_b[L4.y])
local true_persist=_b[L1.y]+_b[L2.y]+_b[L3.y]+_b[L4.y]

xtreg dem i.year l(1/4).y, fe r cluster(wbcode2)
local true_demshock=e(sigma_e)
local true_vardem=e(sigma_u)
local true_dar1=_b[L1.y]
local true_dar2=_b[L2.y]
local true_dar3=_b[L3.y]
local true_dar4=_b[L4.y]

sum dem if year==1960
local true_deminit=r(sd)
sum y if year==1960
local true_gdpinit=r(sd)
/*
/* Simulation assuming initial conditions independent from GDP process*/
simulate ,  saving(${project}/appendix/mcarlo_gmm_nstat_0963, replace) reps(1000) seed(12345): mcsimul, stationary(0) dimT(46) dimN(175) period(4) beta(`true_beta') ar1(`true_ar1') ar2(`true_ar2') ar3(`true_ar3') ar4(`true_ar4') dar1(`true_dar1') dar2(`true_dar2') dar3(`true_dar3') dar4(`true_dar4')  vardem(`true_vardem') vargdp(`true_vargdp') gdpshock(`true_gdpshock') demshock(`true_demshock') gdpinit(`true_gdpinit') deminit(`true_deminit')
gen rel_bias_pers=-1+_b_persistence/`true_persist'
gen rel_bias_short=-1+_b_shortrun/`true_beta' 
gen rel_bias_tstat=-1+_b_shortrun_t/`true_beta_t' 
gen rel_bias_longrun=-1+_b_longrun/`true_longest'
sum _b* rel_bias*

/* Imposes persistence of 0.97*/
local true_ar1=`true_ar1'*0.97/`true_persist'
local true_ar2=`true_ar2'*0.97/`true_persist'
local true_ar3=`true_ar3'*0.97/`true_persist'
local true_ar4=`true_ar4'*0.97/`true_persist'
local true_persist=0.97
local true_longest=`true_beta'/(1-`true_persist')
simulate ,  saving(${project}/appendix/mcarlo_gmm_nstat_097, replace) reps(1000) seed(12345): mcsimul, stationary(0) dimT(46) dimN(175) period(4) beta(`true_beta') ar1(`true_ar1') ar2(`true_ar2') ar3(`true_ar3') ar4(`true_ar4') dar1(`true_dar1') dar2(`true_dar2') dar3(`true_dar3') dar4(`true_dar4')  vardem(`true_vardem') vargdp(`true_vargdp') gdpshock(`true_gdpshock') demshock(`true_demshock') gdpinit(`true_gdpinit') deminit(`true_deminit')
gen rel_bias_pers=-1+_b_persistence/`true_persist'
gen rel_bias_short=-1+_b_shortrun/`true_beta' 
gen rel_bias_tstat=-1+_b_shortrun_t/`true_beta_t' 
gen rel_bias_longrun=-1+_b_longrun/`true_longest'
sum _b* rel_bias*
*/


/* Imposes persistence of 0.98*/
local true_ar1=`true_ar1'*0.98/`true_persist'
local true_ar2=`true_ar2'*0.98/`true_persist'
local true_ar3=`true_ar3'*0.98/`true_persist'
local true_ar4=`true_ar4'*0.98/`true_persist'
local true_persist=0.98
local true_longest=`true_beta'/(1-`true_persist')
simulate ,  saving(${project}/appendix/mcarlo_gmm_nstat_098, replace) reps(1000) seed(12345): mcsimul, stationary(0) dimT(46) dimN(175) period(4) beta(`true_beta') ar1(`true_ar1') ar2(`true_ar2') ar3(`true_ar3') ar4(`true_ar4') dar1(`true_dar1') dar2(`true_dar2') dar3(`true_dar3') dar4(`true_dar4')  vardem(`true_vardem') vargdp(`true_vargdp') gdpshock(`true_gdpshock') demshock(`true_demshock') gdpinit(`true_gdpinit') deminit(`true_deminit')
gen rel_bias_pers=-1+_b_persistence/`true_persist'
gen rel_bias_short=-1+_b_shortrun/`true_beta' 
gen rel_bias_tstat=-1+_b_shortrun_t/`true_beta_t' 
gen rel_bias_longrun=-1+_b_longrun/`true_longest'
sum _b* rel_bias*



/* Imposes persistence of 0.99*/
local true_ar1=`true_ar1'*0.99/`true_persist'
local true_ar2=`true_ar2'*0.99/`true_persist'
local true_ar3=`true_ar3'*0.99/`true_persist'
local true_ar4=`true_ar4'*0.99/`true_persist'
local true_persist=0.99
local true_longest=`true_beta'/(1-`true_persist')
simulate ,  saving(${project}/appendix/mcarlo_gmm_nstat_099, replace) reps(1000) seed(12345): mcsimul, stationary(0) dimT(46) dimN(175) period(4) beta(`true_beta') ar1(`true_ar1') ar2(`true_ar2') ar3(`true_ar3') ar4(`true_ar4') dar1(`true_dar1') dar2(`true_dar2') dar3(`true_dar3') dar4(`true_dar4')  vardem(`true_vardem') vargdp(`true_vargdp') gdpshock(`true_gdpshock') demshock(`true_demshock') gdpinit(`true_gdpinit') deminit(`true_deminit')
gen rel_bias_pers=-1+_b_persistence/`true_persist'
gen rel_bias_short=-1+_b_shortrun/`true_beta' 
gen rel_bias_tstat=-1+_b_shortrun_t/`true_beta_t' 
gen rel_bias_longrun=-1+_b_longrun/`true_longest'
sum _b* rel_bias*


