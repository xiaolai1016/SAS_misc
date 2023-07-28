
%INCLUDE 'E:\5. SAS\Non-parametric\One-way layout\LD_CI.SAS';
options ls =100 ps= 1000 nodate nocenter nonumber;
data a;
input
plot   var	rep	sev	sub;
datalines;

;
run;
/*Check the dataset*/

proc print data=a;
run;

/* A 1-way analysis for severity ratings. */
/*Before using Proc Mixed, one needs the ranks of the observations*/

proc rank data=a out=a;
var sev; *requests ranks for disease ratings;
ranks r; *ranks are stored under the variable r;
run;

/*One-way analysis with Proc Mixed*/

ods rtf file='your file name here .rtf';
proc mixed data=a anovaf;
title1 '1-way analysis using MIXED';
class var;
model r = var / chisq ;
repeated / type=un(1) group=var;
lsmeans var/pdiff;
run;
ods rtf close;
quit;

ods rtf file='your file name here .rtf';
title1 '1-way analysis using macro; each observation is a subject';
%ld_ci (data=a,var=sev,group=var,alpha=0.05,subject=sub);
run;
ods rtf close;
quit;

