options ls=80 ps=60;
title1 'Radial growth ';

data Radial growth;
input sub Trial temp  stage$  Isolate$   Rep audpc;
datalines;

;

proc print;
run;

proc glm;
     class temp stage isolate;
	 model audpc= temp|stage|isolate;
     means  temp|stage|isolate /lsd line;
	 run;

