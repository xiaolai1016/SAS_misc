options pageno=1;
data factor;
input a $ b $ c $ rep	yield;
datalines;

;;
proc print;
title 'Printout of Factorial Data';
run;

proc mixed method=type3;
class a b c;
model yield=a|b|c;
*lsmeans a b c/pdiff=0.05;
*means a*b*c;
title 'ANOVA assuming A and B are fixed effects';
run;