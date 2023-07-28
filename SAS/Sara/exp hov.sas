data sara;
input exp rep trt audpc;
datalines;

;
run;

proc glm data=sara;
class exp;
model audpc=exp;
means exp/ hovtest=levene(type=abs);
run;

