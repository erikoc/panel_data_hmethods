/**********************/
/*Clean and set memory*/
/**********************/
clear all
set matsize 5000
set more off

/******************/
/*Install programs*/
/******************/
cap ssc install xtabond2 
cap ssc install xtivreg2 
cap ssc install spmat 
cap ssc install spmack

/*********************/
/*Sets base directory*/
/*********************/
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */
global limit=25                                                                             /* Evaluate effects 25 years after transition */
local repsBS=100                                                                            /* Number of bootstrap repetitions            */

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/*********************DEFINE REQUIRED PROGRAMS THAT WILL BE USED DURING THE EXECUTION OF THIS DO FILE **************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
capture program drop vareffects
program define vareffects, eclass

quietly: nlcom (effect1: _b[shortrun]) ///
	  (demo1: _b[democonstant]) ///
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post

quietly: nlcom (effect2: _b[effect1]*_b[lag1]+_b[shortrun]) ///
	  (demo2: _b[demo1]*_b[lag1d]+_b[democonstant]) /// 
	  (effect1: _b[effect1]) ///
	  (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post

quietly: nlcom (effect3: _b[effect2]*_b[lag1]+_b[effect1]*_b[lag2]+_b[shortrun]) ///
	  (demo3: _b[demo2]*_b[lag1d]+_b[demo1]*_b[lag2d]+_b[democonstant]) /// 
	  (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post
	  
quietly: nlcom (effect4: _b[effect3]*_b[lag1]+_b[effect2]*_b[lag2]+_b[effect1]*_b[lag3]+_b[shortrun]) ///
	  (demo4: _b[demo3]*_b[lag1d]+_b[demo2]*_b[lag2d]+_b[demo1]*_b[lag3d]+_b[democonstant]) /// 
	  (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post
	  
quietly: nlcom (effect5: _b[effect4]*_b[lag1]+_b[effect3]*_b[lag2]+_b[effect2]*_b[lag3]+_b[effect1]*_b[lag4]+_b[shortrun]) ///
	  (demo5: _b[demo4]*_b[lag1d]+_b[demo3]*_b[lag2d]+_b[demo2]*_b[lag3d]+_b[demo1]*_b[lag4d]+_b[democonstant]) /// 
	  (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post
	  
quietly: nlcom (effect6: _b[effect5]*_b[lag1]+_b[effect4]*_b[lag2]+_b[effect3]*_b[lag3]+_b[effect2]*_b[lag4]+_b[effect1]*_b[lag5]+_b[shortrun]) ///
	  (demo6: _b[demo5]*_b[lag1d]+_b[demo4]*_b[lag2d]+_b[demo3]*_b[lag3d]+_b[demo2]*_b[lag4d]+_b[demo1]*_b[lag5d]+_b[democonstant]) /// 
	  (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect7: _b[effect6]*_b[lag1]+_b[effect5]*_b[lag2]+_b[effect4]*_b[lag3]+_b[effect3]*_b[lag4]+_b[effect2]*_b[lag5]+_b[effect1]*_b[lag6]+_b[shortrun]) ///
	  (demo7: _b[demo6]*_b[lag1d]+_b[demo5]*_b[lag2d]+_b[demo4]*_b[lag3d]+_b[demo3]*_b[lag4d]+_b[demo2]*_b[lag5d]+_b[demo1]*_b[lag6d]+_b[democonstant]) /// 
	  (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect8: _b[effect7]*_b[lag1]+_b[effect6]*_b[lag2]+_b[effect5]*_b[lag3]+_b[effect4]*_b[lag4]+_b[effect3]*_b[lag5]+_b[effect2]*_b[lag6]+_b[effect1]*_b[lag7]+_b[shortrun]) ///
	  (demo8: _b[demo7]*_b[lag1d]+_b[demo6]*_b[lag2d]+_b[demo5]*_b[lag3d]+_b[demo4]*_b[lag4d]+_b[demo3]*_b[lag5d]+_b[demo2]*_b[lag6d]+_b[demo1]*_b[lag7d]+_b[democonstant]) /// 
	  (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	
	  
quietly: nlcom (effect9: _b[effect8]*_b[lag1]+_b[effect7]*_b[lag2]+_b[effect6]*_b[lag3]+_b[effect5]*_b[lag4]+_b[effect4]*_b[lag5]+_b[effect3]*_b[lag6]+_b[effect2]*_b[lag7]+_b[effect1]*_b[lag8]+_b[shortrun]) ///
	  (demo9: _b[demo8]*_b[lag1d]+_b[demo7]*_b[lag2d]+_b[demo6]*_b[lag3d]+_b[demo5]*_b[lag4d]+_b[demo4]*_b[lag5d]+_b[demo3]*_b[lag6d]+_b[demo2]*_b[lag7d]+_b[demo1]*_b[lag8d]+_b[democonstant]) /// 
	  (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	
	  	  
	  
quietly: nlcom (effect10: _b[effect9]*_b[lag1]+_b[effect8]*_b[lag2]+_b[effect7]*_b[lag3]+_b[effect6]*_b[lag4]+_b[effect5]*_b[lag5]+_b[effect4]*_b[lag6]+_b[effect3]*_b[lag7]+_b[effect2]*_b[lag8]+_b[shortrun]) ///
	  (demo10: _b[demo9]*_b[lag1d]+_b[demo8]*_b[lag2d]+_b[demo7]*_b[lag3d]+_b[demo6]*_b[lag4d]+_b[demo5]*_b[lag5d]+_b[demo4]*_b[lag6d]+_b[demo3]*_b[lag7d]+_b[demo2]*_b[lag8d]+_b[democonstant]) /// 
	  (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect11: _b[effect10]*_b[lag1]+_b[effect9]*_b[lag2]+_b[effect8]*_b[lag3]+_b[effect7]*_b[lag4]+_b[effect6]*_b[lag5]+_b[effect5]*_b[lag6]+_b[effect4]*_b[lag7]+_b[effect3]*_b[lag8]+_b[shortrun]) ///
	  (demo11: _b[demo10]*_b[lag1d]+_b[demo9]*_b[lag2d]+_b[demo8]*_b[lag3d]+_b[demo7]*_b[lag4d]+_b[demo6]*_b[lag5d]+_b[demo5]*_b[lag6d]+_b[demo4]*_b[lag7d]+_b[demo3]*_b[lag8d]+_b[democonstant]) /// 
	  (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

	  
quietly: nlcom (effect12: _b[effect11]*_b[lag1]+_b[effect10]*_b[lag2]+_b[effect9]*_b[lag3]+_b[effect8]*_b[lag4]+_b[effect7]*_b[lag5]+_b[effect6]*_b[lag6]+_b[effect5]*_b[lag7]+_b[effect4]*_b[lag8]+_b[shortrun]) ///
	  (demo12: _b[demo11]*_b[lag1d]+_b[demo10]*_b[lag2d]+_b[demo9]*_b[lag3d]+_b[demo8]*_b[lag4d]+_b[demo7]*_b[lag5d]+_b[demo6]*_b[lag6d]+_b[demo5]*_b[lag7d]+_b[demo4]*_b[lag8d]+_b[democonstant]) /// 
	  (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect13: _b[effect12]*_b[lag1]+_b[effect11]*_b[lag2]+_b[effect10]*_b[lag3]+_b[effect9]*_b[lag4]+_b[effect8]*_b[lag5]+_b[effect7]*_b[lag6]+_b[effect6]*_b[lag7]+_b[effect5]*_b[lag8]+_b[shortrun]) ///
	  (demo13: _b[demo12]*_b[lag1d]+_b[demo11]*_b[lag2d]+_b[demo10]*_b[lag3d]+_b[demo9]*_b[lag4d]+_b[demo8]*_b[lag5d]+_b[demo7]*_b[lag6d]+_b[demo6]*_b[lag7d]+_b[demo5]*_b[lag8d]+_b[democonstant]) /// 
	  (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	
	  
quietly: nlcom (effect14: _b[effect13]*_b[lag1]+_b[effect12]*_b[lag2]+_b[effect11]*_b[lag3]+_b[effect10]*_b[lag4]+_b[effect9]*_b[lag5]+_b[effect8]*_b[lag6]+_b[effect7]*_b[lag7]+_b[effect6]*_b[lag8]+_b[shortrun]) ///
	  (demo14: _b[demo13]*_b[lag1d]+_b[demo12]*_b[lag2d]+_b[demo11]*_b[lag3d]+_b[demo10]*_b[lag4d]+_b[demo9]*_b[lag5d]+_b[demo8]*_b[lag6d]+_b[demo7]*_b[lag7d]+_b[demo6]*_b[lag8d]+_b[democonstant]) /// 
	  (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	
	  
quietly: nlcom (effect15: _b[effect14]*_b[lag1]+_b[effect13]*_b[lag2]+_b[effect12]*_b[lag3]+_b[effect11]*_b[lag4]+_b[effect10]*_b[lag5]+_b[effect9]*_b[lag6]+_b[effect8]*_b[lag7]+_b[effect7]*_b[lag8]+_b[shortrun]) ///
	  (demo15: _b[demo14]*_b[lag1d]+_b[demo13]*_b[lag2d]+_b[demo12]*_b[lag3d]+_b[demo11]*_b[lag4d]+_b[demo10]*_b[lag5d]+_b[demo9]*_b[lag6d]+_b[demo8]*_b[lag7d]+_b[demo7]*_b[lag8d]+_b[democonstant]) /// 
	  (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect16: _b[effect15]*_b[lag1]+_b[effect14]*_b[lag2]+_b[effect13]*_b[lag3]+_b[effect12]*_b[lag4]+_b[effect11]*_b[lag5]+_b[effect10]*_b[lag6]+_b[effect9]*_b[lag7]+_b[effect8]*_b[lag8]+_b[shortrun]) ///
	  (demo16: _b[demo15]*_b[lag1d]+_b[demo14]*_b[lag2d]+_b[demo13]*_b[lag3d]+_b[demo12]*_b[lag4d]+_b[demo11]*_b[lag5d]+_b[demo10]*_b[lag6d]+_b[demo9]*_b[lag7d]+_b[demo8]*_b[lag8d]+_b[democonstant]) /// 
	  (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	
	  
quietly: nlcom (effect17: _b[effect16]*_b[lag1]+_b[effect15]*_b[lag2]+_b[effect14]*_b[lag3]+_b[effect13]*_b[lag4]+_b[effect12]*_b[lag5]+_b[effect11]*_b[lag6]+_b[effect10]*_b[lag7]+_b[effect9]*_b[lag8]+_b[shortrun]) ///
	  (demo17: _b[demo16]*_b[lag1d]+_b[demo15]*_b[lag2d]+_b[demo14]*_b[lag3d]+_b[demo13]*_b[lag4d]+_b[demo12]*_b[lag5d]+_b[demo11]*_b[lag6d]+_b[demo10]*_b[lag7d]+_b[demo9]*_b[lag8d]+_b[democonstant]) /// 
	  (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect18: _b[effect17]*_b[lag1]+_b[effect16]*_b[lag2]+_b[effect15]*_b[lag3]+_b[effect14]*_b[lag4]+_b[effect13]*_b[lag5]+_b[effect12]*_b[lag6]+_b[effect11]*_b[lag7]+_b[effect10]*_b[lag8]+_b[shortrun]) ///
	  (demo18: _b[demo17]*_b[lag1d]+_b[demo16]*_b[lag2d]+_b[demo15]*_b[lag3d]+_b[demo14]*_b[lag4d]+_b[demo13]*_b[lag5d]+_b[demo12]*_b[lag6d]+_b[demo11]*_b[lag7d]+_b[demo10]*_b[lag8d]+_b[democonstant]) /// 
	  (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect19: _b[effect18]*_b[lag1]+_b[effect17]*_b[lag2]+_b[effect16]*_b[lag3]+_b[effect15]*_b[lag4]+_b[effect14]*_b[lag5]+_b[effect13]*_b[lag6]+_b[effect12]*_b[lag7]+_b[effect11]*_b[lag8]+_b[shortrun]) ///
	  (demo19: _b[demo18]*_b[lag1d]+_b[demo17]*_b[lag2d]+_b[demo16]*_b[lag3d]+_b[demo15]*_b[lag4d]+_b[demo14]*_b[lag5d]+_b[demo13]*_b[lag6d]+_b[demo12]*_b[lag7d]+_b[demo11]*_b[lag8d]+_b[democonstant]) /// 
	  (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect20: _b[effect19]*_b[lag1]+_b[effect18]*_b[lag2]+_b[effect17]*_b[lag3]+_b[effect16]*_b[lag4]+_b[effect15]*_b[lag5]+_b[effect14]*_b[lag6]+_b[effect13]*_b[lag7]+_b[effect12]*_b[lag8]+_b[shortrun]) ///
	  (demo20: _b[demo19]*_b[lag1d]+_b[demo18]*_b[lag2d]+_b[demo17]*_b[lag3d]+_b[demo16]*_b[lag4d]+_b[demo15]*_b[lag5d]+_b[demo14]*_b[lag6d]+_b[demo13]*_b[lag7d]+_b[demo12]*_b[lag8d]+_b[democonstant]) /// 
	  (effect19: _b[effect19]) (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo19: _b[demo19]) (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect21: _b[effect20]*_b[lag1]+_b[effect19]*_b[lag2]+_b[effect18]*_b[lag3]+_b[effect17]*_b[lag4]+_b[effect16]*_b[lag5]+_b[effect15]*_b[lag6]+_b[effect14]*_b[lag7]+_b[effect13]*_b[lag8]+_b[shortrun]) ///
	  (demo21: _b[demo20]*_b[lag1d]+_b[demo19]*_b[lag2d]+_b[demo18]*_b[lag3d]+_b[demo17]*_b[lag4d]+_b[demo16]*_b[lag5d]+_b[demo15]*_b[lag6d]+_b[demo14]*_b[lag7d]+_b[demo13]*_b[lag8d]+_b[democonstant]) /// 
	  (effect20: _b[effect20]) (effect19: _b[effect19]) (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo20: _b[demo20]) (demo19: _b[demo19]) (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect22: _b[effect21]*_b[lag1]+_b[effect20]*_b[lag2]+_b[effect19]*_b[lag3]+_b[effect18]*_b[lag4]+_b[effect17]*_b[lag5]+_b[effect16]*_b[lag6]+_b[effect15]*_b[lag7]+_b[effect14]*_b[lag8]+_b[shortrun]) ///
	  (demo22: _b[demo21]*_b[lag1d]+_b[demo20]*_b[lag2d]+_b[demo19]*_b[lag3d]+_b[demo18]*_b[lag4d]+_b[demo17]*_b[lag5d]+_b[demo16]*_b[lag6d]+_b[demo15]*_b[lag7d]+_b[demo14]*_b[lag8d]+_b[democonstant]) /// 
	  (effect21: _b[effect21]) (effect20: _b[effect20]) (effect19: _b[effect19]) (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo21: _b[demo21]) (demo20: _b[demo20]) (demo19: _b[demo19]) (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect23: _b[effect22]*_b[lag1]+_b[effect21]*_b[lag2]+_b[effect20]*_b[lag3]+_b[effect19]*_b[lag4]+_b[effect18]*_b[lag5]+_b[effect17]*_b[lag6]+_b[effect16]*_b[lag7]+_b[effect15]*_b[lag8]+_b[shortrun]) ///
	  (demo23: _b[demo22]*_b[lag1d]+_b[demo21]*_b[lag2d]+_b[demo20]*_b[lag3d]+_b[demo19]*_b[lag4d]+_b[demo18]*_b[lag5d]+_b[demo17]*_b[lag6d]+_b[demo16]*_b[lag7d]+_b[demo15]*_b[lag8d]+_b[democonstant]) /// 
	  (effect22: _b[effect22]) (effect21: _b[effect21]) (effect20: _b[effect20]) (effect19: _b[effect19]) (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo22: _b[demo22]) (demo21: _b[demo21]) (demo20: _b[demo20]) (demo19: _b[demo19]) (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect24: _b[effect23]*_b[lag1]+_b[effect22]*_b[lag2]+_b[effect21]*_b[lag3]+_b[effect20]*_b[lag4]+_b[effect19]*_b[lag5]+_b[effect18]*_b[lag6]+_b[effect17]*_b[lag7]+_b[effect16]*_b[lag8]+_b[shortrun]) ///
	  (demo24: _b[demo23]*_b[lag1d]+_b[demo22]*_b[lag2d]+_b[demo21]*_b[lag3d]+_b[demo20]*_b[lag4d]+_b[demo19]*_b[lag5d]+_b[demo18]*_b[lag6d]+_b[demo17]*_b[lag7d]+_b[demo16]*_b[lag8d]+_b[democonstant]) /// 
	  (effect23: _b[effect23]) (effect22: _b[effect22]) (effect21: _b[effect21]) (effect20: _b[effect20]) (effect19: _b[effect19]) (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo23: _b[demo23]) (demo22: _b[demo22]) (demo21: _b[demo21]) (demo20: _b[demo20]) (demo19: _b[demo19]) (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

quietly: nlcom (effect25: _b[effect24]*_b[lag1]+_b[effect23]*_b[lag2]+_b[effect22]*_b[lag3]+_b[effect21]*_b[lag4]+_b[effect20]*_b[lag5]+_b[effect19]*_b[lag6]+_b[effect18]*_b[lag7]+_b[effect17]*_b[lag8]+_b[shortrun]) ///
	  (demo25: _b[demo24]*_b[lag1d]+_b[demo23]*_b[lag2d]+_b[demo22]*_b[lag3d]+_b[demo21]*_b[lag4d]+_b[demo20]*_b[lag5d]+_b[demo19]*_b[lag6d]+_b[demo18]*_b[lag7d]+_b[demo17]*_b[lag8d]+_b[democonstant]) /// 
	  (effect24: _b[effect24]) (effect23: _b[effect23]) (effect22: _b[effect22]) (effect21: _b[effect21]) (effect20: _b[effect20]) (effect19: _b[effect19]) (effect18: _b[effect18]) (effect17: _b[effect17]) (effect16: _b[effect16]) (effect15: _b[effect15]) (effect14: _b[effect14]) (effect13: _b[effect13]) (effect12: _b[effect12]) (effect11: _b[effect11]) (effect10: _b[effect10]) (effect9: _b[effect9]) (effect8: _b[effect8]) (effect7: _b[effect7]) (effect6: _b[effect6]) (effect5: _b[effect5]) (effect4: _b[effect4]) (effect3: _b[effect3]) (effect2: _b[effect2]) (effect1: _b[effect1]) ///
	  (demo24: _b[demo24]) (demo23: _b[demo23]) (demo22: _b[demo22]) (demo21: _b[demo21]) (demo20: _b[demo20]) (demo19: _b[demo19]) (demo18: _b[demo18]) (demo17: _b[demo17]) (demo16: _b[demo16]) (demo15: _b[demo15]) (demo14: _b[demo14]) (demo13: _b[demo13]) (demo12: _b[demo12]) (demo11: _b[demo11]) (demo10: _b[demo10]) (demo9: _b[demo9]) (demo8: _b[demo8]) (demo7: _b[demo7]) (demo6: _b[demo6]) (demo5: _b[demo5]) (demo4: _b[demo4]) (demo3: _b[demo3]) (demo2: _b[demo2]) (demo1: _b[demo1]) /// 
	  (democonstant: _b[democonstant]) (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) (lag2: _b[lag2]) (lag3: _b[lag3]) (lag4: _b[lag4]) (lag5: _b[lag5]) (lag6: _b[lag6]) (lag7: _b[lag7]) (lag8: _b[lag8]) ///
	  (lag1d: _b[lag1d]) (lag2d: _b[lag2d]) (lag3d: _b[lag3d]) (lag4d: _b[lag4d]) (lag5d: _b[lag5d]) (lag6d: _b[lag6d]) (lag7d: _b[lag7d]) (lag8d: _b[lag8d]) ///
	  , post	

	  
quietly: nlcom (effect25: _b[effect25]) (demo25: _b[demo25]) ///
               (counter: _b[effect1]*(_b[demo24]-_b[demo23])+ ///
			   _b[effect2]*(_b[demo23]-_b[demo22])+ ///
			   _b[effect3]*(_b[demo22]-_b[demo21])+ ///
			   _b[effect4]*(_b[demo21]-_b[demo20])+ ///
			   _b[effect5]*(_b[demo20]-_b[demo19])+ ///
			   _b[effect6]*(_b[demo19]-_b[demo18])+ ///
			   _b[effect7]*(_b[demo18]-_b[demo17])+ ///
			   _b[effect8]*(_b[demo17]-_b[demo16])+ ///
			   _b[effect9]*(_b[demo16]-_b[demo15])+ ///
			   _b[effect10]*(_b[demo15]-_b[demo14])+ ///
			   _b[effect11]*(_b[demo14]-_b[demo13])+ ///
			   _b[effect12]*(_b[demo13]-_b[demo12])+ ///
			   _b[effect13]*(_b[demo12]-_b[demo11])+ ///
			   _b[effect14]*(_b[demo11]-_b[demo10])+ ///
			   _b[effect15]*(_b[demo10]-_b[demo9])+ ///
			   _b[effect16]*(_b[demo9]-_b[demo8])+ ///
			   _b[effect17]*(_b[demo8]-_b[demo7])+ ///
			   _b[effect18]*(_b[demo7]-_b[demo6])+ ///
			   _b[effect19]*(_b[demo6]-_b[demo5])+ ///
			   _b[effect20]*(_b[demo5]-_b[demo4])+ ///
			   _b[effect21]*(_b[demo4]-_b[demo3])+ ///
	           _b[effect22]*(_b[demo3]-_b[demo2])+ ///
			   _b[effect23]*(_b[demo2]-_b[demo1])+ ///
			   _b[effect24]*_b[demo1])  ///
			   (effect25counter: _b[effect25]-(_b[effect1]*(_b[demo24]-_b[demo23])+ ///
			   _b[effect2]*(_b[demo23]-_b[demo22])+ ///
			   _b[effect3]*(_b[demo22]-_b[demo21])+ ///
			   _b[effect4]*(_b[demo21]-_b[demo20])+ ///
			   _b[effect5]*(_b[demo20]-_b[demo19])+ ///
			   _b[effect6]*(_b[demo19]-_b[demo18])+ ///
			   _b[effect7]*(_b[demo18]-_b[demo17])+ ///
			   _b[effect8]*(_b[demo17]-_b[demo16])+ ///
			   _b[effect9]*(_b[demo16]-_b[demo15])+ ///
			   _b[effect10]*(_b[demo15]-_b[demo14])+ ///
			   _b[effect11]*(_b[demo14]-_b[demo13])+ ///
			   _b[effect12]*(_b[demo13]-_b[demo12])+ ///
			   _b[effect13]*(_b[demo12]-_b[demo11])+ ///
			   _b[effect14]*(_b[demo11]-_b[demo10])+ ///
			   _b[effect15]*(_b[demo10]-_b[demo9])+ ///
			   _b[effect16]*(_b[demo9]-_b[demo8])+ ///
			   _b[effect17]*(_b[demo8]-_b[demo7])+ ///
			   _b[effect18]*(_b[demo7]-_b[demo6])+ ///
			   _b[effect19]*(_b[demo6]-_b[demo5])+ ///
			   _b[effect20]*(_b[demo5]-_b[demo4])+ ///
			   _b[effect21]*(_b[demo4]-_b[demo3])+ ///
	           _b[effect22]*(_b[demo3]-_b[demo2])+ ///
			   _b[effect23]*(_b[demo2]-_b[demo1])+ ///
			   _b[effect24]*_b[demo1])), post	

			   
	  ereturn display
end

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/************************************ ESTIMATION PROCEDURES AND STORING RESULTS ************************************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/

use ${project}/DDCGdata_final.dta
file open myfile using ${project}/results/TableMain_lags.tex, write replace /*p-value for lags 5 to 8 not reported in the table*/
file write myfile "p-value lags 5 to 8"

/*******************************
Left panel: Within estimator
*******************************/

fvset base 2000 year
fvset base 1 wbcode2

quietly: reg y l.y dem i.year i.wbcode2
estimates store e_y
quietly: reg dem l.dem
estimates store e_dem
quietly: suest e_y e_dem, cluster(wbcode2)
estimates store e1
nlcom (shortrun: _b[e_y_mean:dem])  (lag1: _b[e_y_mean:L.y])  (lag2: 0)  (lag3: 0)  (lag4: 0)   (lag5: 0)   (lag6: 0)   (lag7: 0)   (lag8: 0) ///
      (democonstant: _b[e_dem_mean:_cons])  (lag1d: _b[e_dem_mean:L.dem])  (lag2d: 0)  (lag3d: 0)  (lag4d: 0)   (lag5d: 0)   (lag6d: 0)   (lag7d: 0)   (lag8d: 0), post
vareffects
estimates store e1add

quietly: reg y l(1/2).y dem i.year i.wbcode2
estimates store e_y
quietly: reg dem l(1/2).dem
estimates store e_dem
quietly: suest e_y e_dem, cluster(wbcode2)
estimates store e2
nlcom (shortrun: _b[e_y_mean:dem])  (lag1: _b[e_y_mean:L.y])  (lag2: _b[e_y_mean:L2.y])  (lag3: 0)  (lag4: 0)   (lag5: 0)   (lag6: 0)   (lag7: 0)   (lag8: 0) ///
      (democonstant: _b[e_dem_mean:_cons])  (lag1d: _b[e_dem_mean:L.dem])  (lag2d: _b[e_dem_mean:L2.dem])  (lag3d: 0)  (lag4d: 0)   (lag5d: 0)   (lag6d: 0)   (lag7d: 0)   (lag8d: 0), post
vareffects
estimates store e2add

quietly: reg y l(1/4).y dem i.year i.wbcode2
estimates store e_y
quietly: reg dem l(1/4).dem
estimates store e_dem
quietly: suest e_y e_dem, cluster(wbcode2)
estimates store e3
nlcom (shortrun: _b[e_y_mean:dem])  (lag1: _b[e_y_mean:L.y])  (lag2: _b[e_y_mean:L2.y])  (lag3:  _b[e_y_mean:L3.y])  (lag4:  _b[e_y_mean:L4.y])   (lag5: 0)   (lag6: 0)   (lag7: 0)   (lag8: 0) ///
      (democonstant: _b[e_dem_mean:_cons])  (lag1d: _b[e_dem_mean:L.dem])  (lag2d: _b[e_dem_mean:L2.dem])  (lag3d: _b[e_dem_mean:L3.dem])  (lag4d: _b[e_dem_mean:L4.dem])   (lag5d: 0)   (lag6d: 0)   (lag7d: 0)   (lag8d: 0), post
vareffects
estimates store e3add

quietly: reg y l(1/8).y dem i.year i.wbcode2
estimates store e_y
quietly: reg dem l(1/8).dem
estimates store e_dem
quietly: suest e_y e_dem, cluster(wbcode2)
estimates store e4
nlcom (shortrun: _b[e_y_mean:dem])  (lag1: _b[e_y_mean:L.y])  (lag2: _b[e_y_mean:L2.y])  (lag3:  _b[e_y_mean:L3.y])  (lag4:  _b[e_y_mean:L4.y])   (lag5: _b[e_y_mean:L5.y])   (lag6: _b[e_y_mean:L6.y])   (lag7: _b[e_y_mean:L7.y])   (lag8: _b[e_y_mean:L8.y]) ///
      (democonstant: _b[e_dem_mean:_cons])  (lag1d: _b[e_dem_mean:L.dem])  (lag2d: _b[e_dem_mean:L2.dem])  (lag3d: _b[e_dem_mean:L3.dem])  (lag4d: _b[e_dem_mean:L4.dem])   (lag5d: _b[e_dem_mean:L5.dem])   (lag6d: _b[e_dem_mean:L6.dem])   (lag7d: _b[e_dem_mean:L7.dem])   (lag8d: _b[e_dem_mean:L8.dem]), post
vareffects
estimates store e4add


estout e1 e2 e3 e4  using ${project}/appendix/TableCounter.tex, style(tex) /// 
varlabels(e_y_mean:L.y "log GDP first lag" e_y_mean:L2.y "log GDP second lag" e_y_mean:L3.y "log GDP third lag" e_y_mean:L4.y "log GDP fourth lag" e_y_mean:L5.y "log GDP fifth lag" e_y_mean:L6.y "log GDP sixth lag" e_y_mean:L7.y "log GDP seventh lag" e_y_mean:L8.y "log GDP eight lag"  e_y_mean:dem "Democracy" /// 
e_dem_mean:L.dem "Democracy first lag" e_dem_mean:L2.dem "Democracy second lag" e_dem_mean:L3.dem "Democracy third lag" e_dem_mean:L4.dem "Democracy fourth lag" e_dem_mean:L5.dem "Democracy fifth lag" e_dem_mean:L6.dem "Democracy first lag" e_dem_mean:L7.dem "Democracy seventh lag" e_dem_mean:L8.dem "Democracy eight lag" e_dem_mean:_cons "Constant") ///
cells(b(star fmt(%9.3f)) se(par)) ///
keep(e_y_mean:L.y e_y_mean:L2.y e_y_mean:L3.y e_y_mean:L4.y e_y_mean:L5.y e_y_mean:L6.y e_y_mean:L7.y e_y_mean:L8.y  e_y_mean:dem e_dem_mean:L.dem e_dem_mean:L2.dem e_dem_mean:L3.dem e_dem_mean:L4.dem e_dem_mean:L5.dem e_dem_mean:L6.dem e_dem_mean:L7.dem e_dem_mean:L8.dem e_dem_mean:_cons) ///
order(e_y_mean:dem e_y_mean:L.y e_y_mean:L2.y e_y_mean:L3.y e_y_mean:L4.y e_y_mean:L5.y e_y_mean:L6.y e_y_mean:L7.y e_y_mean:L8.y   e_dem_mean:_cons e_dem_mean:L.dem e_dem_mean:L2.dem e_dem_mean:L3.dem e_dem_mean:L4.dem e_dem_mean:L5.dem e_dem_mean:L6.dem e_dem_mean:L7.dem e_dem_mean:L8.dem) ///
stardrop(e_y_mean:L.y e_y_mean:L2.y e_y_mean:L3.y e_y_mean:L4.y e_y_mean:L5.y e_y_mean:L6.y e_y_mean:L7.y e_y_mean:L8.y  e_y_mean:dem e_dem_mean:L.dem e_dem_mean:L2.dem e_dem_mean:L3.dem e_dem_mean:L4.dem e_dem_mean:L5.dem e_dem_mean:L6.dem e_dem_mean:L7.dem e_dem_mean:L8.dem e_dem_mean:_cons) /// /// 
nolabel replace mlabels(none) collabels(none)



estout e1add e2add e3add e4add  using ${project}/appendix/TableCounter_add.tex, style(tex) /// 
varlabels(effect25 "Effect of democracy after 25 years"  counter "Counterfactual growth" demo25 "Counterfactual likelihood of democracy" effect25counter "Effect of democracy relative to counterfactual") ///
cells(b(star fmt(%9.3f)) se(par)) ///
keep(effect25 demo25 counter effect25counter) /// 
order(effect25 demo25 counter effect25counter) ///
stardrop(effect25 demo25 counter effect25counter) ///
nolabel replace mlabels(none) collabels(none)
