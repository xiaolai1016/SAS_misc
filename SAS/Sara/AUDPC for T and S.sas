
data a;
input tem  ls  iso  rep day dis;
trt= tem*100+ls*10+iso;
datalines;


;;

*segarea is the equation that will calculate the AUDPC;
 
proc sort data=a;
by rep trt;
run;

data b;
set a;
segarea=((dis+lag(dis))/2)*(day-lag(day));
if day=0 then segarea=.;
run;

proc summary nway data=b;
class rep trt;
var dis segarea;
output out=new sum(segarea)=audpc max(dis)=ymax;
proc print;
run;


