options pageno=1;
data RCBD;
input SEEDRATE REP YIELD;
datalines;
1	1	9
1	2	9
1	3	10
1	4	8
2	1	9
2	2	9
2	3	8
2	4	9
3	1	10
3	2	9
3	3	10
3	4	9
4	1	0
4	2	0
4	3	1
4	4	0
5	1	8
5	2	7
5	3	6
5	4	4
6	1	7
6	2	5
6	3	5
6	4	9
1	1	9
1	2	8
1	3	.
1	4	.
2	1	9
2	2	9
2	3	8
2	4	5
3	1	10
3	2	8
3	3	10
3	4	10
4	1	0
4	2	0
4	3	1
4	4	0
5	1	8
5	2	4
5	3	8
5	4	8
6	1	0
6	2	8
6	3	8
6	4	6

		






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