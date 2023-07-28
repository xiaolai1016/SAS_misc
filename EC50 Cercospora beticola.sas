options ls=132 ps=60 mprint symbolgen;

%global ds dim x1 x2 x3 x4 x5 x6;

*** Enter the name of the raw data file below.  This assumes that you ***
*** have exported the data ffrom Excel as a tab-delimited text file   ***
*** which has the extension .txt --- obviously you can replace the    ***
*** data set names in the infile with the actual data file name if    ***
*** you prefer.                                                       ***;
%let ds=example;  *<====== Raw data file name goes here.  ****; 

%let dim=6;  *** Most of your data files have 3 concentrations ...    ***;
             *** If you have more, the dimension needs to be changed. ***;
             
*** The actual concentrations used are assigned via Let statements to ***
*** macro variables below.  If your concentrations are different, you ***
*** need to change these values.                                      ***;
%let x1=-4;  %let x2=-3; 	%let x3=-2;  	%let x4=-1; %let x5=0; %let x6=1;

title1 "ED50 Estimation ---    Tests";
					
******* Linear Splines Interpolation Method Used Below - Fu-Chih Cheng Code;
data lsi;
  infile "&ds..txt" delimiter='09'x missover firstobs=4;
  informat isolate $16.;
  input Isolate Rep growth1-growth6 pt1-pt&dim;
  	 if pt1<50 and pt2>=50 then
        ed50 = (50*(&x1-&x2) - pt2*(&x1-&x2) + (pt1-pt2)*&x2)/(pt1-pt2);
     else if pt2<50 and pt3>=50 then
        ed50 = (50*(&x2-&x3) - pt3*(&x2-&x3) + (pt2-pt3)*&x3)/(pt2-pt3);
     else if pt3<50.0 and pt4>=50.0 then
        ed50 = (50*(&x3-&x4) - pt4*(&x3-&x4) + (pt3-pt4)*&x4)/(pt3-pt4);
	else if pt4<50.0 and pt5>=50.0 then
        ed50 = (50*(&x4-&x5) - pt5*(&x4-&x5) + (pt4-pt5)*&x5)/(pt4-pt5);
	else if pt5<50.0 and pt6>=50.0 then
        ed50 = (50*(&x5-&x6) - pt6*(&x5-&x6) + (pt5-pt6)*&x6)/(pt5-pt6);
	 	 else  ed50 =.;
  ed50back=10**ed50;
run;
  
proc print;
  format ed50back 10.4;
		title3 "Data from &ds";
  run;
