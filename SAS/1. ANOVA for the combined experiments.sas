options formdlim='*';
data pico;
input expt trt	rep	value;
cards;	

;
proc anova;
class expt trt rep;
model value= expt rep(expt) trt trt*expt;
means trt/lsd;
title 'anova for treatment comparison';
run;

proc ttest;
class expt;
var value;
run;

