data raw;
title 'Dante commercial lines data 2009 and 2010';
input exp rep trt audpc;
cards;
1	1	1	10
1	1	2	5
1	1	3	0
1	1	4	5
1	1	5	10
1	1	6	5
1	1	7	5
1	1	8	40
1	1	9	20
1	1	10	40
1	1	11	0
1	1	12	15
1	1	13	30
1	1	14	45
1	1	15	5
1	1	16	20
1	1	17	30
1	1	18	20
1	1	19	0
1	1	20	0
1	2	1	10
1	2	2	20
1	2	3	10
1	2	4	0
1	2	5	5
1	2	6	30
1	2	7	5
1	2	8	15
1	2	9	10
1	2	10	15
1	2	11	20
1	2	12	10
1	2	13	0
1	2	14	15
1	2	15	5
1	2	16	15
1	2	17	25
1	2	18	25
1	2	19	10
1	2	20	0
1	3	1	5
1	3	2	5
1	3	3	10
1	3	4	10
1	3	5	0
1	3	6	10
1	3	7	5
1	3	8	0
1	3	9	0
1	3	10	5
1	3	11	10
1	3	12	15
1	3	13	5
1	3	14	10
1	3	15	10
1	3	16	35
1	3	17	20
1	3	18	20
1	3	19	5
1	3	20	0
1	4	1	15
1	4	2	5
1	4	3	10
1	4	4	5
1	4	5	0
1	4	6	20
1	4	7	10
1	4	8	10
1	4	9	10
1	4	10	10
1	4	11	15
1	4	12	25
1	4	13	10
1	4	14	0
1	4	15	5
1	4	16	20
1	4	17	25
1	4	18	10
1	4	19	0
1	4	20	5
2	5	1	10
2	5	2	5
2	5	3	5
2	5	4	15
2	5	5	0
2	5	6	0
2	5	7	5
2	5	8	5
2	5	9	15
2	5	10	10
2	5	11	10
2	5	12	30
2	5	13	0
2	5	14	15
2	5	15	35
2	5	16	30
2	5	17	45
2	5	18	45
2	5	19	5
2	5	20	0
2	6	1	10
2	6	2	20
2	6	3	0
2	6	4	0
2	6	5	10
2	6	6	20
2	6	7	5
2	6	8	0
2	6	9	15
2	6	10	5
2	6	11	5
2	6	12	5
2	6	13	0
2	6	14	5
2	6	15	30
2	6	16	20
2	6	17	25
2	6	18	25
2	6	19	5
2	6	20	5
2	7	1	15
2	7	2	0
2	7	3	10
2	7	4	5
2	7	5	5
2	7	6	15
2	7	7	5
2	7	8	10
2	7	9	15
2	7	10	5
2	7	11	5
2	7	12	5
2	7	13	5
2	7	14	0
2	7	15	35
2	7	16	30
2	7	17	40
2	7	18	45
2	7	19	0
2	7	20	10
2	8	1	0
2	8	2	15
2	8	3	0
2	8	4	20
2	8	5	5
2	8	6	15
2	8	7	0
2	8	8	20
2	8	9	5
2	8	10	5
2	8	11	5
2	8	12	15
2	8	13	5
2	8	14	5
2	8	15	25
2	8	16	20
2	8	17	45
2	8	18	25
2	8	19	0
2	8	20	5

;;
 
proc sort data=raw;
by exp rep trt;
run;
proc glm data=raw;
      class exp rep trt;
      model audpc = exp rep trt;
      lsmeans trt; *<--these are the means of the ranks, same as in Table above;
Contrast	'R vs S'	trt	1	-1	1	-1	1	-1	1	-1	1	-1	1	-1	1	-1	1	-1	1	-1	1	-1;
Contrast	'5g vs 7g'	trt	1	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
Contrast	'7g vs 14g'	trt	0	0	0	0	1	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0;
Contrast	'5g vs 14g'	trt	1	1	0	0	0	0	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0;
Contrast	'5g vs (5g + foli appl)'	trt	1	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
Contrast	'7g vs(7g + foli appl)'	trt	0	0	0	0	1	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0;
Contrast	'14g vs (14g + foli appl)'	trt	0	0	0	0	0	0	0	0	1	1	-1	-1	0	0	0	0	0	0	0	0;
Contrast	'In-fur vs 5g '	trt	1	1	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0	0	0	0	0;
Contrast	'In-fur vs 7g '	trt	0	0	0	0	1	1	0	0	0	0	0	0	-1	-1	0	0	0	0	0	0;
Contrast	'In-fur vs 14g' 	trt	0	0	0	0	0	0	0	0	1	1	0	0	-1	-1	0	0	0	0	0	0;
Contrast	'In-fur vs foli appl'	trt	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	1	1	0	0	0	0;
Contrast	 'foli appl vs 5 g'	trt	1	1	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0	0	0;
Contrast	 'foli appl vs 7g'	trt	0	0	0	0	1	1	0	0	0	0	0	0	0	0	-1	-1	0	0	0	0;
Contrast	 'foli appl vs 14g'	trt	0	0	0	0	0	0	0	0	1	1	0	0	0	0	-1	-1	0	0	0	0;
Contrast	'foli appl vs (5g + foli appl)'	trt	0	0	1	1	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0	0	0;
Contrast	'foli appl vs(7g + foli appl)'	trt	0	0	0	0	0	0	1	1	0	0	0	0	0	0	-1	-1	0	0	0	0;
Contrast	'foli appl vs (14g + foli appl)'	trt	0	0	0	0	0	0	0	0	0	0	1	1	0	0	-1	-1	0	0	0	0;
Contrast	'Non-Inoc check vs foliar'	trt	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	-1	-1;
Contrast	'Inoc check vs  foli appl'	trt	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	-1	0	0;
Contrast	'Non-Inoc check vs 5g'	trt	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1;
Contrast	'Non-Inoc check vs 7g'	trt	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1;
Contrast	'Non-Inoc check vs 14g'	trt	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	-1	-1;
Contrast	'Inoc check vs 5g '	trt	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0;
Contrast	'Inoc check vs 7g '	trt	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0;
Contrast	'Inoc check vs 14g '	trt	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	-1	-1	0	0;
Contrast	'In-fur vs (5g + foli appl)'	trt	0	0	1	1	0	0	0	0	0	0	0	0	-1	-1	0	0	0	0	0	0;
Contrast	'In-fur vs(7g + foli appl)'	trt	0	0	0	0	0	0	1	1	0	0	0	0	-1	-1	0	0	0	0	0	0;
Contrast	'In-fur vs (14g + foli appl)'	trt	0	0	0	0	0	0	0	0	0	0	1	1	-1	-1	0	0	0	0	0	0;
Contrast	'Non-Inoc check vs (5g + foli appl)'	trt	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1;
Contrast	'Non-Inoc check vs(7g + foli appl)'	trt	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	-1	-1;
Contrast	'Non-Inoc check vs (14g + foli appl)'	trt	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	-1	-1;
Contrast	'Inoc check vs (5g + foli appl)'	trt	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0;
Contrast	'Inoc checkvs(7g + foli appl)'	trt	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	-1	-1	0	0;
Contrast	'Inoc check vs (14g + foli appl)'	trt	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	-1	-1	0	0;
Contrast	'In-fur vs Non-Inoc check'	trt	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	-1	-1;
Contrast	'In-fur vs Inoc-check'	trt	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	-1	-1	0	0;
