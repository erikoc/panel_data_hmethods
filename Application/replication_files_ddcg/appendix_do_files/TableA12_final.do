/**********************/
/*Clean and set memory*/
/**********************/
clear all
set matsize 5000
set more off

/******************/
/*Install programs*/
/******************/
cap ssc install spmat 
cap ssc install spmack

/*********************/
/*Sets base directory*/
/*********************/
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/*********************DEFINE REQUIRED PROGRAMS THAT WILL BE USED DURING THE EXECUTION OF THIS DO FILE **************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/

use ${project}/DDCGdata_final.dta
cap ssc install shp2dta
shp2dta using "${project}/geo_data/shapecountries", database(${project}/geo_data/borders_countries) coordinates(${project}/geo_data/coordinates_countries) replace

forvalues j=1960(1)2010{
preserve
keep if year==`j'
keep if cen_lon!=.&cen_lat!=.&y!=.&dem!=.
replace unrest=0 if unrest==.
keep wbcode2 cen_lon cen_lat y dem year unrest _ID
sort wbcode2

spmat idistance dmat cen_lon cen_lat, id(wbcode2) dfunction(dhaversine) norm(row) replace
spmat lag  y_w dmat y
spmat lag  dem_w dmat dem
spmat lag  unrest_w dmat unrest

spmat contiguity dmat_c using ${project}/geo_data/coordinates_countries, id(_ID)  norm(row) replace
spmat lag  y_wc dmat_c y
spmat lag  dem_wc dmat_c dem
spmat lag  unrest_wc dmat_c unrest

keep y_w dem_w unrest_w y_wc dem_wc unrest_wc  wbcode2 year
tempfile m
save `m', replace
restore

merge 1:1 wbcode2 year using `m', update
drop _m

}

sort wbcode2 year

/*********************************/
/*Jackniffed regional averages   */    
/*********************************/

sort wbcode2 year
gen region2=region+demext

replace unrest=0 if unrest==.

foreach var in dem y unrest {
bysort region2 year: egen total=total(`var')
gen jtotal=total-`var'
replace jtotal=total if `var'==.
bysort region2 year: egen count=count(`var')
gen jcount=count-1
replace jcount=count if `var'==.
gen `var'_reg_jk=jtotal/jcount
drop total count jtotal jcount
}

sort wbcode2 year

/**********************************/
/*Diffusion patterns for democracy*/
/**********************************/

/*Lagged regressions*/
xtreg dem L.dem L.dem_reg_jk yy* if y!=., fe r
estimates store e1

xtreg dem L.dem L.dem_w yy* if y!=., fe r
estimates store e2

xtreg dem L.dem L.dem_wc yy* if y!=., fe r
estimates store e3

xtreg dem L.dem L.dem_reg_jk L.dem_w yy* if y!=., fe r
estimates store e4

xtreg dem L.dem L.dem_reg_jk L.dem_wc yy* if y!=., fe r
estimates store e5

xtreg dem L.dem L.dem_reg_jk L.dem_w L.dem_wc yy* if y!=., fe r
estimates store e6


estout e1 e2 e3 e4 e5 e6 using ${project}/appendix/TableDiffusionDem.tex, ///
style(tex) varlabels(L.dem "Lagged democracy" L.dem_reg_jk "Lagged regional democracy"  L.dem_w "Lagged distance-weighted democracy" L.dem_wc "Lagged neighbors' average democracy") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N N_g, fmt(%7.0f %7.0f) labels( " Observations" "Countries in sample"))   collabels(none) ///
keep(L.dem L.dem_reg_jk L.dem_w L.dem_wc) order(L.dem L.dem_reg_jk L.dem_w L.dem_wc) stardrop(L.dem L.dem_reg_jk L.dem_w L.dem_wc) nolabel replace mlabels(none)


/**********************************/
/*Diffusion patterns for unrest   */
/**********************************/

/*Lagged regressions*/
xtreg unrest L.unrest L.unrest_reg_jk yy* if y!=. , fe r
estimates store e1

xtreg unrest L.unrest L.unrest_w yy* if y!=. , fe r
estimates store e2

xtreg unrest L.unrest L.unrest_wc yy* if y!=. , fe r
estimates store e3

xtreg unrest L.unrest L.unrest_reg_jk L.unrest_w yy* if y!=. , fe r
estimates store e4

xtreg unrest L.unrest L.unrest_reg_jk L.unrest_wc yy* if y!=. , fe r
estimates store e5

xtreg unrest L.unrest L.unrest_reg_jk L.unrest_w L.unrest_wc yy* if y!=. , fe r
estimates store e6

estout e1 e2 e3 e4 e5 e6 using ${project}/appendix/TableDiffusionUnrest.tex, ///
style(tex) varlabels(L.unrest "Lagged unrest" L.unrest_reg_jk "Lagged regional unrest"  L.unrest_w "Lagged distance-weighted unrest" L.unrest_wc "Lagged neighbors' average unrest") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N N_g, fmt(%7.0f %7.0f) labels( " Observations" "Countries in sample"))   collabels(none) ///
keep(L.unrest L.unrest_reg_jk L.unrest_w L.unrest_wc) order(L.unrest L.unrest_reg_jk L.unrest_w L.unrest_wc) stardrop(L.unrest L.unrest_reg_jk L.unrest_w L.unrest_wc) nolabel replace mlabels(none)


/**********************************/
/*Diffusion patterns for Income   */
/**********************************/

/*Lagged regressions*/
xtreg y L.y L.y_reg_jk yy* if y!=. , fe r
estimates store e1

xtreg y L.y L.y_w yy* if y!=. , fe r
estimates store e2

xtreg y L.y L.y_wc yy* if y!=. , fe r
estimates store e3

xtreg y L.y L.y_reg_jk L.y_w yy* if y!=. , fe r
estimates store e4

xtreg y L.y L.y_reg_jk L.y_wc yy* if y!=. , fe r
estimates store e5

xtreg y L.y L.y_reg_jk L.y_w L.y_wc yy* if y!=. , fe r
estimates store e6

estout e1 e2 e3 e4 e5 e6  using ${project}/appendix/TableDiffusionGDP.tex, ///
style(tex) varlabels(L.y "Lagged GDP" L.y_reg_jk "Lagged regional GDP"  L.y_w "Lagged distance-weighted GDP" L.y_wc "Lagged neighbors' average GDP") ///
cells(b(star fmt(%9.3f)) se(par)) stats( N N_g, fmt(%7.0f %7.0f) labels( " Observations" "Countries in sample"))   collabels(none) ///
keep(L.y L.y_reg_jk L.y_w L.y_wc) order(L.y L.y_reg_jk L.y_w L.y_wc) stardrop(L.y L.y_reg_jk L.y_w L.y_wc) nolabel replace mlabels(none)



