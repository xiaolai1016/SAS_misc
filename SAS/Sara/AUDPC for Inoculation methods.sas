
data a;
input iso method rep day dis;
trt=iso*10+method;
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


