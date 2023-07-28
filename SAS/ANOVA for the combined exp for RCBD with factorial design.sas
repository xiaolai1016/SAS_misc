options formdlim='*';
data pico;
input expt a b $	rep	value;
cards;	
options pageno=1;
data RCBD;
input SEEDRATE REP YIELD;
datalines;

;;

proc print;
title 'Printout of RCBD Data with No Sampling';
run;

proc anova;
classes rep seedrate;
model yield=rep seedrate;
means seedrate/lsd;
title 'ANOVA for RCBD with No Sampling';
run;


proc anova;
class expt a b rep ;
model value= expt|a|b rep(expt);
test h=a e=expt*a;
test h=b e=expt*b;
test h=a*b e=expt*a*b;
means a /lsd e=expt*a;
means b /lsd e=expt*b;
means a*b/lsd e=expt*a*b lines;
title 'anova for treatment comparison';
run;

proc ttest;
class expt;
var value;
run;

proc means stderr;
class a b;
var value;
run;