
%INCLUDE 'C:\Users\Xiao\Desktop\Data\5. SAS\Non-parametric\One-way layout\LD_CI.SAS';
options ls =100 ps= 1000 nodate nocenter nonumber;
data sara;

input rep   a b sev ;
treatment=10*a+b;
subject=100*rep+treatment;
datalines;

;；

run;

proc print data=sara;
run;



proc rank data=sara out=sara;
var sev; *requests ranks for disease ratings;
ranks r; *ranks are stored under the variable r;
run;


proc sort data=sara out=sara;
by a b; 
run;

proc mixed data=sara anovaf method=MIVQUE0;
class a b;
model r=a b a*b/chisq;
repeated/type=un(1) group=a*b;
lsmeans a b a*b;
run; 

%ld_ci (data=sara,var=sev,group=treatment,alpha=0.05,subject=subject);
run;

ods rtf close;
quit;

