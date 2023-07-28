
options formdlim='*';
data Homogenity;
input expt	trt	rep	severity;
cards; 

;

proc anova;
class expt trt rep;
model severity=expt;
means expt/hovtest=BARTLETT;
title 'Barellets test';
run;
