%include 'C:\Yangxi\Data\5. SAS\Non-parametric\One-way layout\LD_CI.SAS';
%INCLUDE 'C:\Yangxi\Data\5. SAS\Non-parametric\Split-plot\F1_LD_F1.SAS';

data sara;
input rep tem ls iso  dis;
trt=100*tem + 10*ls + iso;
sub=1000*rep + trt;
datalines;

;
run;

proc rank data=sara out=sara;
var dis;
ranks r;
run;

proc sort data=sara out=sara;
by tem ls iso;
run;

proc mixed data=sara;
class tem ls iso;
model r=tem|ls|iso/chisq;
repeated/type=un(1) group=tem*ls*iso;
lsmeans tem|ls|iso;
run;
%LD_CI(data=sara, var=dis, group=trt, alpha=0.05, subject=sub);
run; 

