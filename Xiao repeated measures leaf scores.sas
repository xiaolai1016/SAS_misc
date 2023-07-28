/*
Experiment conducted by Xiao Lai on methods of incoulation of Fusarium on sugar beet
run= two runs (experiment conducted twice)
method = dipping drench cutting barroot barseed
isol = Fus_oxys and Fus_sec
rep = 4 replications 
dai= days after inoculation
s = disease severity rating 0-5 scale on 0 7 14 21 28 and 35 days after inoculation
*/
dm 'log;clear;output;clear;';
%INCLUDE 'F:\My Vaults\literature\Statistical analyses\Non-Parametric analysis\APS Workshop Non-parametric statistical analysis\Macros\LD_CI.sas';
%INCLUDE 'F:\My Vaults\literature\Statistical analyses\Non-Parametric analysis\APS Workshop Non-parametric statistical analysis\Macros\F2_LD_F1.sas';
%INCLUDE 'F:\My Vaults\literature\Statistical analyses\Non-Parametric analysis\APS Workshop Non-parametric statistical analysis\Macros\mult.sas';

options ls =100 ps= 1000 nodate nonumber;
data xiao;
input run isol method rep dai s0 s7 s14 s21 s28 s35;
group = right(method)||"_"||left(put(isol,2.)); *<--Combination of variety and isol;
sub1 = _n_; *<--Unique label for every inoc x fert x rep combination;
cards;
1	1	1	1	0	0	0	1	5	5	5
1	1	1	2	0	0	0	0	5	5	5
1	1	1	3	0	0	0	5	5	5	5
1	1	1	4	0	0	0	0	5	5	5
1	1	1	5	0	0	0	3	4	5	5
1	1	1	6	0	0	0	0	0	5	5
1	1	2	1	0	0	0	5	0	0	0
1	1	2	2	0	0	0	5	0	0	0
1	1	2	3	0	0	0	5	0	0	0
1	1	2	4	0	0	0	5	1	1	1
1	1	2	5	0	0	0	4	0	0	0
1	1	2	6	0	0	0	0	0	4	5
1	1	3	1	0	0	0	5	1	1	1
1	1	3	2	0	0	0	5	2	5	5
1	1	3	3	0	0	0	5	4	5	5
1	1	3	4	0	0	0	5	1	2	3
1	1	3	5	0	0	0	5	0	2	4
1	1	3	6	0	0	0	5	1	4	5
1	1	4	1	0	0	0	0	0	0	5
1	1	4	2	0	0	0	0	3	3	5
1	1	4	3	0	0	0	0	4	5	5
1	1	4	4	0	0	0	0	3	5	5
1	1	4	5	0	0	0	1	4	5	5
1	1	4	6	0	0	0	0	5	5	5
1	1	5	1	0	0	0	3	5	5	5
1	1	5	2	0	0	0	4	5	5	5
1	1	5	3	0	0	0	3	5	5	5
1	1	5	4	0	0	0	4	5	5	5
1	1	5	5	0	0	0	2	3	4	4
1	1	5	6	0	0	0	4	5	5	5
1	1	6	1	0	0	0	5	5	5	5
1	1	6	2	0	0	0	0	1	4	5
1	1	6	3	0	0	0	2	3	3	3
1	1	6	4	0	0	0	4	5	5	5
1	1	6	5	0	0	0	0	1	4	5
1	1	6	6	0	0	0	0	1	4	5
1	2	1	1	0	0	0	0	3	5	5
1	2	1	2	0	0	0	1	5	5	5
1	2	1	3	0	0	0	4	5	5	5
1	2	1	4	0	0	0	5	5	5	5
1	2	1	5	0	0	0	4	5	5	5
1	2	1	6	0	0	0	0	0	4	5
1	2	2	1	0	0	0	0	0	0	0
1	2	2	2	0	0	0	0	0	0	0
1	2	2	3	0	0	0	0	0	0	0
1	2	2	4	0	0	0	0	0	0	0
1	2	2	5	0	0	0	0	1	1	1
1	2	2	6	0	0	0	0	0	0	0
1	2	3	1	0	0	0	0	0	1	1
1	2	3	2	0	0	0	0	4	4	4
1	2	3	3	0	0	0	0	3	3	4
1	2	3	4	0	0	0	1	3	4	4
1	2	3	5	0	0	0	0	0	2	3
1	2	3	6	0	0	0	0	0	0	1
1	2	4	1	0	0	0	0	3	3	3
1	2	4	2	0	0	0	0	3	3	3
1	2	4	3	0	0	0	0	1	3	3
1	2	4	4	0	0	0	1	3	3	3
1	2	4	5	0	0	0	0	4	4	5
1	2	4	6	0	0	0	0	1	1	1
1	2	5	1	0	0	0	2	5	5	5
1	2	5	2	0	0	0	3	4	5	5
1	2	5	3	0	0	0	3	5	5	5
1	2	5	4	0	0	0	3	4	5	5
1	2	5	5	0	0	0	3	4	4	4
1	2	5	6	0	0	0	2	4	4	5
1	2	6	1	0	0	0	0	0	1	3
1	2	6	2	0	0	0	0	0	3	4
1	2	6	3	0	0	0	0	0	0	3
1	2	6	4	0	0	0	0	3	3	3
1	2	6	5	0	0	0	0	0	3	3
1	2	6	6	0	0	0	0	0	1	3
2	1	1	1	0	0	0	0	3	5	5
2	1	1	2	0	0	0	0	5	5	5
2	1	1	3	0	0	0	0	5	5	5
2	1	1	4	0	0	0	1	5	5	5
2	1	1	5	0	0	0	1	5	5	5
2	1	1	6	0	0	0	5	5	5	5
2	1	2	1	0	0	0	0	0	1	1
2	1	2	2	0	0	0	0	5	5	5
2	1	2	3	0	0	0	0	0	0	1
2	1	2	4	0	0	0	0	0	4	5
2	1	2	5	0	0	0	5	5	5	5
2	1	2	6	0	0	0	0	0	0	1
2	1	3	1	0	0	0	0	0	0	1
2	1	3	2	0	0	0	1	5	5	5
2	1	3	3	0	0	0	0	0	1	3
2	1	3	4	0	0	0	0	1	5	5
2	1	3	5	0	0	0	0	0	0	3
2	1	3	6	0	0	0	0	0	0	3
2	1	4	1	0	0	0	0	0	0	1
2	1	4	2	0	0	0	0	0	0	1
2	1	4	3	0	0	0	0	4	5	5
2	1	4	4	0	0	0	0	5	5	5
2	1	4	5	0	0	0	0	5	5	5
2	1	4	6	0	0	0	0	0	0	0
2	1	5	1	0	0	0	4	5	5	5
2	1	5	2	0	0	0	3	5	5	5
2	1	5	3	0	0	0	4	5	5	5
2	1	5	4	0	0	0	4	5	5	5
2	1	5	5	0	0	0	5	5	5	5
2	1	5	6	0	0	0	2	4	5	5
2	1	6	1	0	0	0	4	5	5	5
2	1	6	2	0	0	0	2	5	5	5
2	1	6	3	0	0	0	4	5	5	5
2	1	6	4	0	0	0	4	5	5	5
2	1	6	5	0	0	0	1	5	5	5
2	1	6	6	0	0	0	3	5	5	5
2	2	1	1	0	0	0	5	5	5	5
2	2	1	2	0	0	0	0	0	4	5
2	2	1	3	0	0	0	3	5	5	5
2	2	1	4	0	0	0	0	0	4	5
2	2	1	5	0	0	0	3	4	5	5
2	2	1	6	0	0	0	0	0	0	2
2	2	2	1	0	0	0	0	0	0	0
2	2	2	2	0	0	0	0	0	0	0
2	2	2	3	0	0	0	0	0	0	0
2	2	2	4	0	0	0	0	0	0	0
2	2	2	5	0	0	0	0	1	3	3
2	2	2	6	0	0	0	0	0	0	0
2	2	3	1	0	0	0	0	0	1	3
2	2	3	2	0	0	0	0	0	1	4
2	2	3	3	0	0	0	0	1	3	3
2	2	3	4	0	0	0	0	0	0	0
2	2	3	5	0	0	0	0	0	0	0
2	2	3	6	0	0	0	0	3	3	3
2	2	4	1	0	0	0	0	0	0	1
2	2	4	2	0	0	0	0	1	1	1
2	2	4	3	0	0	0	0	0	0	0
2	2	4	4	0	0	0	0	0	0	0
2	2	4	5	0	0	0	0	0	3	4
2	2	4	6	0	0	0	2	4	4	4
2	2	5	1	0	0	0	3	5	5	5
2	2	5	2	0	0	0	2	4	5	5
2	2	5	3	0	0	0	2	4	4	4
2	2	5	4	0	0	0	4	5	5	5
2	2	5	5	0	0	0	4	5	5	5
2	2	5	6	0	0	0	3	5	5	5
2	2	6	1	0	0	0	0	1	1	3
2	2	6	2	0	0	0	0	1	1	3
2	2	6	3	0	0	0	0	0	0	1
2	2	6	4	0	0	0	0	1	0	3
2	2	6	5	0	0	0	0	0	1	3
2	2	6	6	0	0	0	0	0	0	3
;
run;

/*First, get the data in univariate format for analysis by macros and Proc Mixed*/
data xiao2;
  set xiao;
  t1=s0;
  t2=s7;
  t3=s14;
  t4=s21;
  t5=s28;
  t6=s35;
  array t{6} t1-t6; *<--Create an array;
  do i=1 to 6;
    if i=1 then time=0;
	if i=2 then time=7;
	if i=3 then time=14;
	if i=4 then time=21;
    if i=5 then time=28;
	if i=6 then time=35;
	sev=t{i};      *<--sev is the variable storing the rating values;
	output;
  end;
  drop i t1-t6 s0 s7 s14 s21 s28 s35; 
run;
title 'Data in univariate format';
proc print data=xiao2;
run;

/*Summary of the median ratings*/
title 'Median ratings';
proc tabulate data=xiao2 f=6.1;
  class method isol time;
  var sev;
  table method*isol, time*sev*(median); *<--Request median rating;
run;

title 'Non-parametric analysis with F2_LD_F1 macro';
%F2_LD_F1(data=xiao2, var=sev, factor1=method, factor2=isol, time=time, subject=sub1);
title 'Relative treatment effects';
%ld_ci(data=xiao2, var=sev, group=group, time=time, subject=sub1);
title 'Marginal effects for methods of inoculation';
%ld_ci(data=xiao2, var=sev, group=method, time=time, subject=sub1);
title 'Marginal effects for isolates';
%ld_ci(data=xiao2, var=sev, group=isol, time=time, subject=sub1);


* Non-parametric analysis usign proc mixed;
/*First, obtain the ranks*/
proc rank data=xiao2 out=xiao3;
var sev;
ranks r;
run;
/*Sort by main factors, then time factor*/
proc sort data=xiao3 out=xiao4;
by sub1 time;
run;
proc mixed data=xiao4 anovaf method=mivque0;
class method isol time sub1;
model r = method|isol|time /chisq;
repeated / type=un sub=sub1 group=method*isol; *<--or group=group, same definition;
lsmeans method*isol*time;
ods output lsmeans=c; *<--Output the lsmeans to dataset c;
run;

* Plotting relative treatment effects. Calculate relative treatment effects from lsmeans;
data c2;
  set c;
  re = (Estimate-0.5)/864; *<-- re are the relative treatment effects. The denominator is the number of observations;
  drop Effect DF tValue Probt; *<-- Just cleaning up a bit;
 proc tabulate data=c2 f=6.2;
  class method isol time;
  var re;
  table isol*method, time*re*(mean); *<--Table of the mean relative effects;
  title 'Estimated treatment relative effects calculated from lsmeans';
run;

/*Now plot re*/
goptions reset=all;
symbol1 v=dot i=join;
title 'RELATIVE TREATMENT EFFECTS';
proc gplot data = c2;
  plot re*time = method;
  plot re*time = isol;
run;
quit;

*To examine in more detail the significant interactions between isolates by time.
One can look at the isol*time relative effects;
title 'EXAMINE THE ISOL*TIME INTERACTION';
proc mixed data=xiao4 anovaf method=mivque0;
class method isol time sub1;
model r = method|isol|time /chisq;
repeated / type=un sub=sub1 group=method*isol;
lsmeans isol*time;
ods output lsmeans=d;
run;

/*                 PLOTTING RELATIVE TREATMENT EFFECTS                */
data d2;
  set d;
  re = (Estimate-0.5)/864; *<-- re are the relative treatment effects;
  drop Effect DF tValue Probt; *<-- Just cleaning up a bit;
run;
title 'RELATIVE TREATMENT EFFECTS FOR ISOL*TIME';
proc tabulate data=d2 f=6.2;
  class isol time;
  var re;
  table isol, time*re*(mean);
run;

/*Now plot re*/
goptions reset=all;
symbol1 v=dot i=join;
title 'RELATIVE TREATMENT EFFECTS BY ISOLATE';
proc gplot data = d2;
  plot re*time = isol;
run;

*To examine in more detail the significant interactions between method by time;
proc mixed data=xiao4 anovaf method=mivque0;
class method isol time sub1;
model r = method|isol|time /chisq;
repeated / type=un sub=sub1 group=method*isol;
lsmeans method*time;
ods output lsmeans=d3;
run;
data d4;
  set d3;
  re = (Estimate-0.5)/864; *<-- re are the relative treatment effects;
  drop Effect DF tValue Probt; *<-- Just cleaning up a bit;
run;
title 'Relative treatment effects for method*time';
proc tabulate data=d4 f=6.2;
  class method time;
  var re;
  table method, time*re*(mean);
run;

/*Now plot re*/
goptions reset=all;
symbol1 v=dot i=join;
title 'Relative treatment effects by method';
proc gplot data = d4;
  plot re*time = method;
run;
quit;


/*One may be interested in all pairwise comparisons between isolates.
This can be done by specifying contrasts.
An alternative is to use the mult macro, written by H.-P. Piepho,
which will produce letter-labelled groupings familiar in Proc GLM output.

NOTE:  use this as a rough estimate only!!
Standard errors produced by the pdiff option in Proc Mixed are not the same as the true standard
errors for the nonparametric analysis.  They may be close enough though to 
use the following as an approximation.  Can double-check with contrast statements.*/
title 'MULT COMPARISONS ON ISOL*TIME';
proc mixed data=xiao4 anovaf method=mivque0;
  class method isol time sub1;
  model r = method|isol|time /chisq;
  repeated / type=un sub=sub1 group=method*isol;
  lsmeans isol*time/pdiff;
  ods output diffs=diffs lsmeans=lsmeans;
%mult(trt=isol, by=time, level=7);
%mult(trt=isol, by=time, level=14);
%mult(trt=isol, by=time, level=21);
%mult(trt=isol, by=time, level=28);
%mult(trt=isol, by=time, level=35);
run;









