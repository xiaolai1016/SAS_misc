options formdlim='*';
data pico;
input expt	trt	rep	value;
cards;	
1	1	1	7
1	1	2	10
1	1	3	7
1	1	4	10
1	2	1	6
1	2	2	9
1	2	3	5
1	2	4	6
1	3	1	10
1	3	2	9
1	3	3	6
1	3	4	10
1	4	1	7
1	4	2	7
1	4	3	6
1	4	4	8
1	5	1	9
1	5	2	9
1	5	3	9
1	5	4	8
1	6	1	0
1	6	2	0
1	6	3	1
1	6	4	0
3	1	1	8
3	1	2	8
3	1	3	7
3	1	4	7
3	2	1	10
3	2	2	10
3	2	3	7
3	2	4	9
3	3	1	8
3	3	2	8
3	3	3	7
3	3	4	9
3	4	1	8
3	4	2	9
3	4	3	9
3	4	4	7
3	5	1	8
3	5	2	9
3	5	3	10
3	5	4	9
3	6	1	1
3	6	2	0
3	6	3	0
3	6	4	0
4	1	1	10
4	1	2	7
4	1	3	8
4	1	4	9
4	2	1	8
4	2	2	9
4	2	3	8
4	2	4	6
4	3	1	10
4	3	2	6
4	3	3	10
4	3	4	9
4	4	1	7
4	4	2	8
4	4	3	9
4	4	4	8
4	5	1	10
4	5	2	9
4	5	3	8
4	5	4	8
4	6	1	2
4	6	2	0
4	6	3	0
4	6	4	0

;
proc anova;
class expt trt rep;
model value= expt| trt ;
means trt*expt/lsd;
title 'anova for treatment comparison';
run;

proc ttest;
class expt;
var value;
run;
