/**********************************************************************/
/*                                                                    */
/*                     F1_LD_F1 Version 1.1                           */
/*                                                                    */
/*                    Februar 2005 fuer SAS Version 9.1               */
/*                                                                    */
/*                    Split - Plot - Design                           */
/*                                                                    */
/* A ist ein fester Faktor mit i=1,...,r unabhaengigen Stufen         */
/* B ist ein fester Faktor mit j=1,...,c   abhaengigen Stufen         */
/* C ist ein zufaelliger Faktor mit k=1,...,n_i unabhaengigen Stufen  */
/*   in jeder Zelle (i,j,k) sind Messwiederholungen erlaubt           */
/* Name: F1_LD_F1.sas                                                 */
/* Quelle : Akritas/Brunner/Langer & Akritas/Brunner                  */
/*          & Brunner/Munzel/Puri                                     */
/* geschrieben von/ written by Steffen Ballerstedt                    */
/* Datum: Januar 1998                                                 */
/* letzte Version: Carola Werner, 1. Februar 2005                     */
/* Änderung, weil SAS 9 bei # (Hadamard) etwas geändert hat.          */ 
/**********************************************************************/

/**********************************************************************/
/* Makroname: contrast                                                */
/* Beschreibung: erzeugt eine Folge von Variablennamen, naemlich      */
/*               C_1,T_1,C_2,T_2,C_3,T_3,C_4,T_4,C_5,T_5              */
/* Aufruf: %contrast                                                  */
/**********************************************************************/

/**********************************************************************/
%macro contrast;
 C_1,T_1
 %do i = 2 %to 5;
  ,C_&i,T_&i
 %end;
%mend;                             /* Ende von contrast               */
/**********************************************************************/

/**********************************************************************/
%macro stepname;
 text1,text2,text3,name1,name2,name3,char1,char2,char3
%mend;
/**********************************************************************/

/**********************************************************************/
%macro F1_LD_F1(
        DATA = , VAR    =    , FACTOR =    , TIME    =,SUBJECT =,
        DATA_PST =_no_, VAR_PST=_no_, FAC_PST=_no_, TIME_PST=_no_,
        DATA_PIT =_no_, VAR_PIT=_no_,               TIME_PIT=_no_,
        DATA_PGT =_no_, VAR_PGT=_no_, FAC_PGT=_no_);
/**********************************************************************/
%LET FACTORA = &FACTOR;
%LET FACTORB = &TIME;
%LET FACTORC = &SUBJECT;
/**********************************************************************/
/* at first it is nessessary to produce a new data-set with the means */
/* of each cell in the variable RmIJK_ = RIJK */
/* sort the data alphanumerical */
proc sort data=&DATA;
 by &FACTORA &FACTORB &FACTORC;

/* calculate the average-ranks in the variable RIJKL */
/* and write it in the data-set daten */
proc rank data=&DATA ties=mean out=daten;
 var &VAR;
 ranks RIJKL;

/* these procedures count the number of levels within factor A*/
/* the data are in the variable count */
proc freq data=daten;
 tables &FACTORA*&FACTORC / out=_hallo_ noprint;
proc freq data=_hallo_;
 tables &FACTORA /out=ni noprint;

/* calculate the means within each cell */
/* missing values are considerd */
/* the data are in the variable  */
proc means data=daten noprint;
 by &FACTORA &FACTORB &FACTORC;
 var RIJKL;
 output out=_value_ N=mIJK MEAN=RmIJK_;

/* clean the data-set middle */
/* this SAS data set contains the means, levels of the factors */
/* and the number of repeated measurement within each cell */
data _value_;
 set _value_;
 keep mIJK RmIJK_ &FACTORB &FACTORA &FACTORC;
run;

/***************************** IML ************************************/
proc iml worksize=120;
reset linesize=80;
/**********************************************************************/
start lese_dat(%stepname,RIJK)global(nI,mIJK);
 use _value_;
  read all var{&FACTORA} into text1;
  read all var{&FACTORB} into text2;
  read all var{&FACTORC} into text3;
  read all var{mIJK}     into mIJK;   /* mIJK is the vector of the    */
                                /* number of measurement in each cell */
  read all var{RmIJK_}   into RIJK;   /* ranks within each cell       */
 close _value_;

 use ni;
  read all var{count} into nI;
 close ni;

finish;
/**********************************************************************/
start set_var(werte,%stepname,index)
       global(r,c,ni,nIJ,mv,mvN,nsum,cnsum,mIJK,N,minj_nij);

 N     = sum(mIJK);                   /* number of total observations */
 cnsum = nrow(text1);                 /* nsum*c=total number of cells */

 x_1   = j(cnsum,1,0);  x_2   = j(cnsum,1,0); x_3 = j(cnsum,1,0);
 name1 = unique(text1); name2 = unique(text2);
 r     = ncol(name1);   c     = ncol(name2);

 /* x_1 berechnen */
 do i=1 to r;
  ind = loc ( text1 = name1[i] );
  x_1[ind] = i;
 end;

 /* x_2 berechnen */
 do i=1 to c;
  ind = loc ( text2 = name2[i] );
  x_2[ind] = i;
 end;

 /* nsum is the total number of subjects */
 nsum = sum(ni);

 /* attention: if one cell is empty (missing data) */
 /* you have to set a . in one row of the data set */
 if (cnsum)^=nrow(mIJK) then print
        'Error',
        'If you have missing values (empty cell(s))',
        'please fill in one observation with a .';

 /* calculate x_3 */
 do i=1 to r;
  m1 = shape((1:ni[i]),c*ni[i],1);
  if i=1 then x_3=m1; else; x_3 = x_3 // m1;
 end;  /* do i */
 /* end of calculation of x_3 */

 index = x_1 || x_2 || x_3;

 run nu_char(name1,char1);
 run nu_char(name2,char2);

 run gen_mv(werte);                  /* generiere die Hilfsmatrix mv  */
 mvN = mv[+,2];                      /* Gesamtanzahl besetzter Zellen */

 nIJ = j(r,c,0);
 do i=1 to r;
  ind     = loc(index[,1]=i*j(cnsum,1,1));
  kinMiJ  = (shape(mv[ind,2],c,ni[i]))`;
  nIJ[i,] = kinMiJ[+,];
  MiJ     = shape(mv[ind,2],c,ni[i]);
  mini    = min(nIJ[i,]);            /* Mindestanzahl der subj ohne   */
/*if i=1 then minj_nij = mini; else;*/
  minj_nij  = minj_nij // mini;          /* missings in jeder Stufe von A */
  /* hier noch keine Kontrolle, ob niju immer > 1 ist              **** ACHTUNG */
 end; /* do i */

finish;
/**********************************************************************/

/**********************************************************************/
start nu_char(vec,char_new);
 if type(vec)='N' then
  do; char_new = left(char(vec)); char_new = right(trim(char_new)); end;
 else char_new = trim(vec);
finish;
/**********************************************************************/

/**********************************************************************/
/* Name: kontrast                                                     */
/* Input-Variablen: keine (global: r,c)                               */
/* Output-Variablen: die Sequenz der Kontrastmatrizen %contrast       */
/* Beschreibung: Berechnung der Kontraste fuer die Hypothesen         */
/*   H_0(A), H_0(B), H_0(AB), H_0(A|B) und H_0(B|A)                   */
/* Nebeneffekte: keine                                                */
/**********************************************************************/

/**********************************************************************/
start kontrast(%contrast)
      global(r,c);

 p_r    = I(r)-j(r,r,1/r) ;
 p_c    = I(c)-j(c,c,1/c) ;

 C_1    = p_r     @  j(1,c,1/c);
 C_2    = j(1,r,1/r)  @  p_c   ;
 C_3    = p_r     @  p_c   ;
 C_4    = p_r     @  I(c)  ;
 C_5    = I(r)    @  p_c   ;

 t_1    = c_1` * ginv(c_1 * c_1`) * c_1;
 t_2    = c_2` * ginv(c_2 * c_2`) * c_2;
 t_3    = c_3` * ginv(c_3 * c_3`) * c_3;
 t_4    = c_4` * ginv(c_4 * c_4`) * c_4;
 t_5    = c_5` * ginv(c_5 * c_5`) * c_5;

finish;                            /* Ende von kontrast               */
/**********************************************************************/

/**********************************************************************/
start gen_mv(werte) global(cnsum,mv);
  mv = j(cnsum,2,1);
  do i=1 to cnsum;
    if werte[i] = . then do;
      mv[i,1] = . ; mv[i,2] = 0 ;
    end;                                  /* do                       */
  end;                                    /* do i                     */
finish;                                   /* Ende von erz_mv          */
/**********************************************************************/

/**********************************************************************/
/* Name: mvspcov                                                      */
/* Input-Variablen: index,RIJK (global: r,c,ni,nIJ,nsum,cnsum,mv,mvN) */
/* Output-Variablen: v_d,RmIJ_                                        */
/* Beschreibung: Berechnung der Kovarianzmatrix gem. Br/Pu/Mu         */
/* Nebeneffekte: keine                                                */
/**********************************************************************/

/**********************************************************************/
start mvspcov(index,RIJK,v_d,RmIJ_)
      global(r,c,ni,nIJ,nsum,cnsum,mv,mvN,N);

 /* v_d mit missings mit RmIJ_ ---------------------------------------*/

 do i=1 to r;                      /* getrennt fuer jede Stufe von A  */
  ind     = loc(index[,1]=i*j(cnsum,1,1));
  RiJKmat = shape(RIJK[ind],c,ni[i]);
  RmiJ_2  = RiJKmat*j(ni[i],1,1)/(nIJ[i,])`;  /* zeilenweise Division */
  if i=1 then RmIJ_ = RmiJ_2; else RmIJ_ = RmIJ_ // RmiJ_2;
  M2      = (shape(RmiJ_2,ni[i],c))`;
  MiJ     = shape(mv[ind,2],c,ni[i]);
  M3      = (RiJKmat - M2) # MiJ;
  niUmat  = shape(nIJ[i,],c,c);
  niJU    = MiJ * MiJ`;
  faktori = 1/((niUmat-1)`#(niUmat-1)+(niJU-1));     /* nach BR/MU/PU */
  *old    = niJU # (1/(niJU-1)) # (1/(niUmat#niUmat`));
  vdi     = M3 * M3` # faktori;    /* Summe ueber k                   */
  if i=1 then v_d = vdi; else v_d = block(v_d,vdi);
 end;                              /* do i                            */

 v_d = v_d * (nsum / (N * N));

finish;                            /* Ende von mvspcov                */
/**********************************************************************/

/**********************************************************************/
start wald_n(anz,p_dach,kontr,v_d,w,fg);
  cp  = kontr * p_dach;
  cvc = kontr * v_d * kontr`;
  w   = anz * cp` * ginv(cvc) * cp;
  fg  = trace(cvc * ginv(cvc));
finish;
/**********************************************************************/

/**********************************************************************/
start box1(p_dach,T,v_d,B1,fg)global(nsum);
 tv = T*v_d;
 B1 = (nsum*p_dach`*T*p_dach)/trace(tv);
 fg = ((trace(tv))**2)/(trace(tv*tv));
finish;  /* box1 */
/**********************************************************************/

/*fuer die simple-time-effects*/
/**********************************************************************/
start box3(ni,p_dach,T,v_d,B1,fg);
 tv = T*v_d;
 B1 = (ni*p_dach`*T*p_dach)/trace(tv);
 fg = ((trace(tv))**2)/(trace(tv*tv));
finish;  /* box3 */
/**********************************************************************/

/**********************************************************************/
start box2(p_dach,T,v_d,B2,fg1,fg2)global(r,c,nsum,ni,cnsum);
 d_p_r = diag(I(r) - j(r,r,1/r));
 mat   = I(r) @ j(1,c,1/c);
 v_a   = mat * v_d *mat`;
 tv    = T*v_d ;
 lambda= inv(diag(ni) - I(r));
 B2    = nsum*(p_dach`*T*p_dach)/(trace(tv));
 fg1   = trace(tv) ** 2 / trace(tv**2);
 fg2   = trace(d_p_r * v_a) **2 / trace(d_p_r**2 * v_a ** 2 * lambda);
finish;

/**********************************************************************/
start pattern(anz,w,p_dach,kontr,v_d,T);
  wc  = w` * kontr;
  cov = sqrt(wc * v_d * wc`);
  T   = (sqrt(anz) * wc * p_dach) / cov;
  /* T ist asymptotisch Standardnormalverteilt */
finish;  /* pattern */
/**********************************************************************/

/***** berechnet die patternd group x time interaction ****************/
start ww_pat(%stepname,p_dach,v_d,tabelle,paare)
      global(c,r,ni,nsum,weights,nIJ);

 sort &DATA_PIT by &TIME_PIT;
 use &DATA_PIT;
  read all var{&VAR_PIT} into weights;
  read all var{&TIME_PIT} into name_pit;
 close &DATA_PIT;

 p_c    = I(c)-j(c,c,1/c) ;
 do i=1 to r-1;
  do j=i+1 to r;
   pos_A    = j(1,r,0); pos_A[i]=1; pos_A[j]=-1;
   con_AB   = pos_A @ p_c;
   run pattern(nsum,weights,p_dach,con_AB,v_d,out_AB);
   pwert_NV = 1 - probnorm(out_AB);

   /* a estimation for the df of the t-distribution */
   mat      = diag(pos_A) @ (weights`*p_c);
   s_n      = mat * v_d * mat`;
   lambda1  = diag(1/(ni-1));
   df       = (trace(s_n))**2/trace(s_n#s_n#lambda1);
   pwert_t  = 1 - probt(out_AB,df);

   nv_and_t = out_AB || pwert_NV || df || pwert_t;
   tabelle  = tabelle // nv_and_t;
   paare    = paare // (char1[i]+'*'+char1[j]);
  end;
 end;

finish; /* ww_pat */
/**********************************************************************/

/**********************************************************************/
/* Name: pst (=pattern-simple-tests)                                  */
/* Input: %stepname (nur: name1), p_dach, v_d                         */
/* Output: name_pa,result2,sim_patm                                   */
/**********************************************************************/

start pst(%stepname,p_dach,v_d,name_pa,result2,sim_patm)
   global(r,c,ni,nsum,cnsum);

 /* Einlesen der pattern ---------------------------------------------*/
 sort &DATA_PST by &FAC_PST &TIME_PST;
 use &DATA_PST;                              /* Einlesen der Pattern  */
  read all var{&VAR_PST}  into sim_pat;      /* fuer die simple tests */
  read all var{&FAC_PST} into pa_text;       /* fuer die simple tests */
 close &DATA_PST;

 /* Analyse, welche der pattern simple tests durchgefuehrt werden ----*/
 name_pa  = unique(pa_text);
 nsimtest = ncol(name_pa);         /* Anzahl der simple-pattern-tests */
 simtest  = j(r,1,0);     /* simtest enthaelt dort 1, wo ein pattern- */
 do z=1 to nsimtest;         /* simple-test durchgefuehrt werden soll */
   ind = loc(name1 = j(1,r,name_pa[z]));
   simtest[ind] = 1;         /* in der Stufe 'ind' des Faktors A soll */
 end;                 /* ein pattern-simple-test durchgefuehrt werden */

 /* die Schleife durch die gesamten Stufen von A ---------------------*/
 do i=1 to r;
             /* das i braucht man fuer die Erzeugung der korrigierten */
                      /* Kovarianzmatrix und des speziellen Kontrasts */

  /* Auswertung der simple-pattern-tests -----------------------------*/
  lauf = 0;        /* Variable, die die Nummer der Auswertung angibt, */
                        /* da unter Umstaenden nicht fuer alle Stufen */
                     /* Gewichte eingegeben wurden und somit auch nur */
               /* die entsprechenden Tests ausgefuehrt werden sollen. */
  if simtest[i] = 1 then do;   /* nur wenn simple test fuer die i-te  */
                               /* Stufe des Faktors A gewuenscht ist  */


   /* die korrigierte Kovarianzmatrix ist simple_v -------------------*/
   simple_v = ni[i] * v_d / nsum;
   einheit= j(r,r,0);
   einheit[i,i] = 1;
   p_c      = I(c)-j(c,c,1/c);                           /* Projektor */
   K_B      = einheit @ p_c;                        /* Kontrastmatrix */
   lauf     = lauf +1;
   sim_patm = shape(sim_pat,sum(simtest),c);
   gewichte = (sim_patm[lauf,])`;      /* der Gewichtsvektor fuer die */
                                                /* lauf-te Auswertung */
   gew      = j(r,1,1) @ gewichte;                      /* aufblaehen */

   /* pattern-Statistik mit Normalverteilungsapproximation -----------*/
   run pattern(ni[i],gew,p_dach,K_B,simple_v,pat_B);
   p_npat   = 1-probnorm(pat_B);                            /* p-Wert */
   npat_B   = pat_B || {.} || p_npat;

   /* pattern-Statistik mit t-Verteilungsapproximation ---------------*/
   fg_tpat  = ni[i] - 1;
   p_tpat   = 1-probt(pat_B,fg_tpat);                       /* p-Wert */
   tpat_B   = pat_B || fg_tpat || p_tpat ;

   /* Untereinanderschreiben der Ergebnisse --------------------------*/
   result2 = result2 // npat_B // tpat_B;

  end;  /* do */
 end;  /* do i */
;

finish;
/**********************************************************************/

/********************* pgt=Pattern-group-test *************************/
start pgt(%stepname,v_d,p_dach,w_pgt,out)global(ni,r,c,nsum);

 out = j(1,4,0);
 sort &DATA_PGT by &FAC_PGT;

 use &DATA_PGT;
  read all var{&VAR_PGT}  into w_pgt;
  read all var{&FAC_PGT}  into name_pgt;
 close &DATA_PGT;

 /*
 %if &DATA_PGT^=_no_ then; do;
 if ((name_pgt`) = name1) then;
  else %put WARNING -> different names in &DATA_PGT and &DATA
 for the levels of the factor &FACTOR;  end;
 */

 /* calc. for the asymptotically normal distributed test statistic */
 p_r    = I(r)-j(r,r,1/r);
 mat    = p_r*(I(r)@j(1,c,1/c));
 run pattern(nsum,w_pgt,p_dach,mat,v_d,T);
 out[1] = T;
 out[2] = 1 - probnorm(T);

 /* estimation for df of the t-distribution */
 ci2    = (diag(w_pgt` * p_r))**2;
 lambda1= diag(1/(ni-1));
 mittel = (I(r)@j(1,c,1/c));
 s_n    = mittel*v_d*mittel`;
 out[3] = (trace(ci2#s_n))**2/(trace((ci2#s_n)**2#lambda1));
 out[4] = 1 - probt(T,out[3]);

finish; /* pgt */
/**********************************************************************/

/*** Funktion fuer die Paarvergleiche innerhalb des Faktors A *********/
start paircomp(%stepname,r,c,v_d,p_dach,paare,test,tabelle);
 p_c = I(c)-j(c,c,1/c) ;
 do i=1 to r-1;
  do j=i+1 to r;
   pos_A    = j(1,r,0); pos_A[i] = 1; pos_A[j] = -1;
   pos_B    = j(1,r,0); pos_B[i] = 1; pos_B[j] = 1;
   con_A    = pos_A   @ j(1,c,1/c);
   con_B    = pos_B/2 @ p_c;
   con_AB   = pos_A   @ p_c;
   t_A      = con_A` * ginv(con_A * con_A`) * con_A;
   t_B      = con_B` * ginv(con_B * con_B`) * con_B;
   t_AB     = con_AB`* ginv(con_AB* con_AB`)* con_AB;
   run box1(p_dach,T_A ,v_d,out_A ,fg_A );
   run box1(p_dach,T_B ,v_d,out_B ,fg_B );
   run box1(p_dach,T_AB,v_d,out_AB,fg_AB);
   if out_A >= 0 & fg_A > 0 then pwert_A  = 1 - probchi(out_A # fg_A, fg_A );
   else pwert_A = .;
   if out_B >= 0 & fg_B > 0 then pwert_B  = 1 - probchi(out_B # fg_B, fg_B );
   else pwert_B = .;
   if out_AB >= 0 & fg_AB > 0 then pwert_AB  = 1 - probchi(out_AB # fg_AB, fg_AB );
   else pwert_AB = .;
   out_A    = out_A || fg_A || pwert_A;
   out_B    = out_B || fg_B || pwert_B;
   out_AB   = out_AB|| fg_AB|| pwert_AB;
   ergebnis = out_A // out_B // out_AB;
   tabelle  = tabelle // ergebnis;
   paare    = paare// j(3,1,(char1[i]+'*'+char1[j])) ;
   /* erzeuge den Vektor mit den Levelnamen --------------------------*/
   if i=1 & j=2
    then test="&FACTORA"//"&FACTORB"//("&FACTORA"+'*'+"&FACTORB");
    else test=test//"&FACTORA"//"&FACTORB"//("&FACTORA"+'*'+"&FACTORB");
  end;
 end;
finish; /* paircomp */
/**********************************************************************/

/****** simple time effekt ********************************************/

start simple_B(i,p_dach,v_d,T1)                                /* Achtung !!!! */
      global(r,c,ni,nsum);                                     /* siehe sim96! */
/* diese simple-tests werden immer fuer alle Stufen von A ausgefuehrt */

/* die korrigierte Kovarianzmatrix ist simple_v */
simple_v = ni[i] * v_d / nsum;

 einheit= j(r,r,0);
 einheit[i,i] = 1;
 p_c    = I(c)-j(c,c,1/c);                               /* Projektor */
 K_B    = einheit @ p_c;                            /* Kontrastmatrix */

/* Wald mit Chi-Quadrat-Approximation --------------------------------*/
 run wald_n(ni[i],p_dach,K_B,simple_v,W_chi,fg_B);
 if W_chi >= 0 & fg_B > 0 then p_w_chi = 1-probchi(W_chi,fg_B);  /* p-Wert */
 else p_w_chi=.;
 w_chi2  = w_chi || fg_B || p_w_chi;

/* ANOVA-Typ-Statistik -----------------------------------------------*/
 run box3(ni[i],p_dach,K_B,simple_v,B1,fg);
 if B1 >= 0 & fg > 0 then p_box1 = 1 - probchi(B1*fg,fg);
 else p_box1 = .;
 w_box1 = B1 || fg || p_box1;

/* Untereinanderschreiben der beiden Ergebnisse ----------------------*/
 T1 = w_chi2 // w_box1;
finish;
/**********************************************************************/

/**********************************************************************/
start rang_vec(RmIJ_,rg_vec)
      global(r,c);

RmIJ_mat = shape(RmIJ_,r,c);
Rm_J_    = RmIJ_mat[+,]/r;                                /* fuer out */
RmI__    = RmIJ_mat[,+]/c;                                /* fuer out */
rg_vec   = RmI__ // Rm_J_` // RmIJ_;                      /* fuer out */

finish;  /* Ende von rang_vec */
/**********************************************************************/

/**********************************************************************/
start tests(RIJK,index,%contrast,%stepname,wald,box,box1,_2x2,
            pat_AB,pat_pair,simple1,RmIJ_,name_pa,result2,sim_patm,v_d,
            paare,pair_tst,pair_tab,tab_pgt,w_pgt)
     global(r,c,nsum,mvN,mv,N,mIJK);

 INDICATOR = loc(mv[,2]=0);
 IF NCOL(INDICATOR) >0 THEN RIJK [INDICATOR]=0; /* set the missings 0,Aenderung 2005*/

 run mvspcov(index,RIJK,v_d,RmIJ_);/* Kovarianz schaetzen             */
 p_dach = (RmIJ_-0.5)/N;           /* Vektor der rel. Beh.effekte     */

 if r=2 & c=2 & mIJK=1 then do;
   R1K1 = RIJK[loc(index[,1] = 1 & index[,2] = 1)]; R11Q = R1K1[:];
   R1K2 = RIJK[loc(index[,1] = 1 & index[,2] = 2)]; R12Q = R1K2[:];
   R2K1 = RIJK[loc(index[,1] = 2 & index[,2] = 1)]; R21Q = R2K1[:];
   R2K2 = RIJK[loc(index[,1] = 2 & index[,2] = 2)]; R22Q = R2K2[:];
   ni = nrow(R1K1) // nrow(R2K1);
   sigmaqn = ssq(R1K1 + R1K2 - R11Q - R12Q) / (ni[1] - 1) / ni[1] //
             ssq(R2K1 + R2K2 - R21Q - R22Q) / (ni[2] - 1) / ni[2];
   tauqn = ssq(R1K1 - R1K2 - R11Q + R12Q) / (ni[1] - 1) / ni[1] //
           ssq(R2K1 - R2K2 - R21Q + R22Q) / (ni[2] - 1) / ni[2];
   una = (R11Q + R12Q - R21Q - R22Q) / sqrt(sigmaqn[+]);
   unt = (R11Q - R12Q + R21Q - R22Q) / sqrt(tauqn[+]);
   unat = (R11Q - R12Q - R21Q + R22Q) / sqrt(tauqn[+]);
   nuea = sigmaqn[+] ** 2 / (sigmaqn ## 2 # (ni - 1) ## -1)[+];
   nuet = tauqn[+] ** 2 / (tauqn ## 2 # (ni - 1) ## -1)[+];
   _2x2 = (una || 2 * probnorm(-abs(una)) || 2 * probt(-abs(una),nuea) || nuea) //
          (unt || 2 * probnorm(-abs(unt)) || 2 * probt(-abs(unt),nuet) || nuet) //
          (unat || 2 * probnorm(-abs(unat)) || 2 * probt(-abs(unat),nuet) || nuet);
 end;
 else do;
    /* Deklarationen ----------------------------------------------------*/
    W_wald = j(3,1,0); fg_wald = j(3,1,0); p_wald  = j(3,1,0);
    B_box1 = j(3,1,0); fg_box1 = j(3,1,0); p_box1  = j(3,1,0);

    /* Auswerten der Statistiken ----------------------------------------*/

    %macro test_mac;
     %do i = 1 %to 3;
      /* wald-type-statistics berechnen ---------------------------------*/
      run wald_n(nsum,p_dach,C_&i,v_d,W,fg);   /* Wald-Typ-Statistik v_d */
      W_wald[&i]  = w;                         /* Statistiken mit v_d    */
      fg_wald[&i] = fg;                        /* allg. Freiheitsgrade   */

      /* anova-type-statistics with box-approximation -------------------*/
      run box1(p_dach,T_&i,v_d, B1,fg);        /* Anova-Typ-Stat mit v_d */
      B_box1[&i]  =  B1;                   /* Statistiken mit v_d    */

   *if fg<0.0001 then kennung=1;                               /**** ACHTUNG ****/
   *if fg<0.0001 then goto neu1;                               /**** ACHTUNG ****/

      fg_box1[&i] = fg;                        /* Freiheitsgrade         */
     %end;
    %mend;                                     /* Ende von test_mac      */
    %test_mac;                                 /* Aufruf                 */

    if W_wald[1] >= 0 & fg_wald[1] > 0 
      then p_wald[1] = 1 - probchi(W_wald[1],fg_wald[1]);
    else p_Wald[2] = .;
    if W_wald[2] >= 0 & fg_wald[2] > 0 
      then p_wald[2] = 1 - probchi(W_wald[2],fg_wald[2]);
    else p_Wald[2] = .;
    if W_wald[3] >= 0 & fg_wald[3] > 0 
      then p_wald[3] = 1 - probchi(W_wald[3],fg_wald[3]);
    else p_Wald[3] = .;

    if B_box1[1] >= 0 & fg_box1[1] > 0 
      then p_box1[1] = 1 - probchi((B_box1#fg_box1)[1],fg_box1[1]);
    else p_box1[1] = .;
    if B_box1[2] >= 0 & fg_box1[2] > 0 
      then p_box1[2] = 1 - probchi((B_box1#fg_box1)[2],fg_box1[2]);
    else p_box1[2] = .;
    if B_box1[3] >= 0 & fg_box1[3] > 0 
      then p_box1[3] = 1 - probchi((B_box1#fg_box1)[3],fg_box1[3]);
    else p_box1[3] = .;

    wald  = W_wald  || fg_wald || p_wald;
    box1  = B_box1  || fg_box1 || p_box1;

    /* anova-type-statistics with double-box-approximation - for H_0(A) -*/
    run box2(p_dach,T_1,v_d,B_box,fg1,fg2);
    if B_box >= 0 & fg1 > 0 & fg2 > 0 then p_box = 1 - probf(B_box,fg1,fg2);
    else p_box = .;
    box   = B_box || fg1 || fg2 || p_box;

    /* simple tests werden fuer alle Stufen von A durchgefuehrt ---------*/
    do i=1 to r;
     run simple_B(i,p_dach,v_d,result1);
     /* i wird fuer die richtige Initialisierung uebergeben */
     if i=1 then simple1 = result1; else simple1 = simple1 // result1;
    end; /* do i */
 end;

 /* der pattern-interaction-test ---------- pairwise -----------------*/
 %if &VAR_PIT ^=_no_ %then %do;
  run ww_pat(%stepname,p_dach,v_d,pat_AB,pat_pair); %end;

 /* der pattern-group-test--------------------------------------------*/
 %if &VAR_PGT^=_no_ %then %do;
  run pgt(%stepname,v_d,p_dach,w_pgt,tab_pgt); %end;

 /* pattern-simple-test(s) nur wenn &DATA_PST vorhanden ---------------*/
 %if &DATA_PST ^= _no_ %then %do;
  %put the Datafile for the pattern-simple-tests is &DATA_PST ;
  run pst(%stepname,p_dach,v_d,name_pa,result2,sim_patm);
  %end; /* %do */

 /* pairwise comparisons fuer mehr als zwei Stufen im Faktor A -------*/
 if r>2
  then run paircomp(%stepname,r,c,v_d,p_dach,paare,pair_tst,pair_tab);

finish; /* tests */
/**********************************************************************/

/**********************************************************************/

start ausgabe(v_d,source_w,source_b,sourceb1,source_2,pat_AB,pat_pair,simple1,
              %stepname,RmIJ_,name_pa,result2,sim_patm,
              paare,pair_tst,pair_tab,tab_pgt,w_pgt)
      global(r,c,ni,nIJ,weights,cnsum,mvN,nsum,N,mIJK);

/*--- Ausgabe allgemeiner Informationen zum Design -------------------*/

 print 'F1_LD_F1 --- subjects(A) x T',
       'A(=FACTOR), T(=TIME): fixed, subjects: random';

/*--- Ausgabe allgemeiner Informationen zum Datensatz ----------------*/

 print "SAS-datafile-name: &data",
       "Response variable: &var";
 print 'Class Level Information' ;
 class  = {A &factorA,T &factorB};
 levels = r // c ;
 print class levels;

 reset noname;
 missing = cnsum - mvN;
 print 'Total number of observations' N,
       '    Number of missing values' missing;

/*--- Output der Rangmittel und der Zellbesetzungen ------------------*/

 %let FAK_A = %trim(&FACTORA) ;            /* Leerzeichen abschneiden */
 %let FAK_B = %trim(&FACTORB) ;            /* Leerzeichen abschneiden */

 spalte1 = j(r,1,"&FAK_A") // j(c,1,"&FAK_B")
           // j(r*c,1,"&FAK_A"+'*'+"&FAK_B"+' ') ;
 do i=1 to r ;
    do j=1 to c ;
       text = text // (char1[i] + '*' + char2[j] )  ;
    end;
 end;
 spalte2 = (char1`) // (char2`)  // text ;
 source = concat (spalte1,spalte2);         /* 2 Spalten verschmelzen */

 run rang_vec(RmIJ_,rg_vec);
 rel_eff = (rg_vec-0.5)/N;
 nor = nIJ[,+] // (nIJ[+,])` // shape(nIJ,r*c,1);

 print 'RTE  = Relative Treatment Effects          ',
       'Nobs = Number of observations (do not count',
       'the repeated measurements within the cells)';
 print source[c={source}] rg_vec[format=6.5 c={"Rank mean"}]
       nor[c={"Nobs"}] rel_eff[c={"RTE"}];

/*--- Ausgabe eventueller Warnungen --------------------------------*/

 if c > nsum | det(v_d) < 10**(-10) then
   print / '';
 reset nocenter;
 if c > nsum then
   print 'Warning:',
         'There are less subjects than sub-plot factor levels.';
 if det(v_d)< 10**(-10) then
   print 'Warning:',
         'Do not use the Wald-type-statistic, because the covariance matrix is singular.';
 *if min(eigval(v_d)) < 0 then
 *  print 'Warning:',
 *        'The estimated covariance matrix is not positive semidefinite.';
 reset center;

/*--- Ausgabe der Testergebnisse -------------------------------------*/
 print /;
 if r=2 & c=2 & mIJK=1 then do;
   zeilen = {'A', 'T', 'AT'};
   spalten   = {'', 'N(0,1)', 't_n', 'n'} ;
   print 'Tests for Group Effect, Time Effect and Interaction (2x2-Design)';
   print'Approximation for large sample sizes with normal-distribution',
        'Approximation for small sample sizes with t_DF';
   print '                           STATISTIC      P_VALUE        DF                    ',
       source_2[format=6.5 r=zeilen c=spalten];
 end;
 else do;
   zeilen = {'A', 'T', 'AT'};
   sp_w   = {W DF P_VALUE} ;
   print 'Wald-type statistic','Approximation for
   large sample sizes with Chi-Square_DF     ';
   print source_w[format=6.5 r=zeilen c=sp_w];

   sp_b1  = {B DF P_VALUE} ;
   print 'ANOVA-type statistic',
   'Box-approximation for small sample sizes with Chi-square_DF';
   print sourceb1[format=6.5 r=zeilen c=sp_b1];

   sp_b   = {B DF1 DF2 P_VALUE} ;
   print 'ANOVA-type statistic',
   'modified Box-approximation for the whole-plot factor A',
   'for small sample sizes with F(DF1,DF2)';
   print source_b[format=6.5 r={A} c=sp_b];
   print /;

  /*--- Ausgabe der Testergebnisse der simple-tests --------------------*/
   sp_simp = {'T' 'DF1' 'P_VALUE'};
   names   = shape({'Wald' 'ANOVA'},2*r,1);
   stufen1 = shape((shape(char1,2,r))`,2*r,1);
   print 'Tests for the simple >>' "&FACTORB" '<< effect (T)',
   'Wald-type (Chi-square_DF1, asymptotic)',
   'ANOVA-type (Chi-square_DF1/DF1, asymptotic)';
   print names[c='Statistic'] stufen1[c="&FAK_A"] simple1[format=6.5 c=sp_simp];
 end;
 print / ;

/*--- Ausgabe der Testergebnisse der pattern-simple-tests ------------*/
 %if &DATA_PST ^= _no_ %then %do;
  run ausgabe2(%stepname,name_pa,result2,sim_patm); %end;

/*--- Ausgabe der Testergebnisse der pairwise comparisons ------------*/
 if r>2 then run ausgabe3(paare,pair_tst,pair_tab);

/*--- Ausgabe der Testergebnisse der pattern interaction tests -------*/
 %if &VAR_PIT^=_no_ %then %do;
  run ausgabe4(%stepname,pat_AB,pat_pair); %end;

/*--- Ausgabe des Testergebnisses des pattern group tests ------------*/
 %if &VAR_PGT^=_no_ %then %do;
  run ausgabe5(%stepname,tab_pgt,w_pgt); %end;

finish;  /* ausgabe */
/**********************************************************************/

/**********************************************************************/
start ausgabe2(%stepname,name_pa,result2,sim_patm);

/*--- Ausgabe der Testergebnisse der pattern-simple-tests ------------*/
*  if type(name_pa) = 'N'
*   then char_pa = char( name_pa,max( int(log10(abs(name_pa+1))) ) );
*   else char_pa = trim(name_pa);
char_pa = name_pa;

 print "Pattern-Statistics (normal,t_DF) for the simple time effects",
       "SAS-Datafile: &DATA_PST";
 print'Approximation for large sample sizes with normal-distribution',
       'Approximation for small sample sizes with t_DF';
   sp_simp2= {'T' 'DF' 'P_VALUE'};
   nsimtest = ncol(name_pa);
   stat2   = shape({'normal' 't'},2*nsimtest,1);
   stufen2 = shape((shape(char_pa,2,nsimtest))`,2*nsimtest,1);
   print stufen2 result2[format=6.5 r=stat2 c=sp_simp2];
   print 'Each row is one pattern';
   print (char_pa`)[c='Levels'] sim_patm[c=char2];
   print / ;
finish; /* ausgabe2 */

/**********************************************************************/
start ausgabe3(paare, test, tabelle);
 print 'Test for pairwise comparisons';
 reset noname;
 print paare[c={pairs}] test[c={test}] tabelle[format=6.5 c={F DF P_VALUE}];
 reset name;
 print /;
finish; /* ausgabe3 */
/**********************************************************************/

/**********************************************************************/
start ausgabe4(%stepname,tabelle,pat_pair) global(weights);

 print"Pattern-Test for pairwise &FactorA * &FactorB interaction";
 print'Approximation for large sample sizes with normal-distribution',
      'Approximation for small sample sizes with t_DF',
      "SAS-Datafile: &DATA_PIT, Pattern-variable: &VAR_PIT";

 reset noname; print pat_pair tabelle[format=6.5 c={T P_VALUE_NV DF P_VALUE_t_DF}];

 levels = weights`;
 reset name; print levels[c=char2 r={pattern}]; reset noname;

finish; /* ausgabe4 */
/**********************************************************************/

/**********************************************************************/

start ausgabe5(%stepname,tabelle,w_pgt);
 print "Pattern-group-test (group=&FACTORA)";
 print'Approximation for large sample sizes with normal-distribution',
      'Approximation for small sample sizes with t_DF',
      "SAS-Datafile: &DATA_PGT, Pattern-variable: &VAR_PGT";

 reset noname; print tabelle[format=6.5 c={T P_VALUE_NV DF P_VALUE_t_DF}];
 levels = w_pgt`;
 reset name; print levels[c=char1 r={pattern}]; reset noname;

finish;

/**********************************************************************/

/**********************************************************************/
run lese_dat(%stepname,RIJK);
run set_var(RIJK,%stepname,index);
run kontrast(%contrast);
run tests(RIJK,index,%contrast,%stepname,wald,
          box,box1,_2x2,pat_AB,pat_pair,simple1,RmIJ_,name_pa,result2,sim_patm,v_d,
          paare,pair_tst,pair_tab,tab_pgt,w_pgt);
run ausgabe(v_d,wald,box,box1,_2x2,pat_AB,pat_pair,simple1,
            %stepname,RmIJ_,name_pa,result2,sim_patm,
            paare,pair_tst,pair_tab,tab_pgt,w_pgt);
/**********************************************************************/


/**********************************************************************/
quit; /* IML */
/**********************************************************************/
%mend F1_LD_F1;
/**********************************************************************/

/* Die Reihenfolge der Stufen innerhalb der Faktoren wird automatisch
   alphanumerisch vorgenommen. */

/*
Der Funktionsablauf:
--------------------
f1_ld_f1  -> lese_dat
          -> set_var -> gen_mv
          -> kontrast
          -> tests   -> mvranks
                     -> mvspcov
                     -> if r=2 and c=2 and all mijk=1:-
                     -> else: -> %test_mac
                                 -> wald_n (5 x)
                                 -> box1   (5 x)
                              -> box2
                              -> simple_B (r x)
                     -> if &VAR_PIT^=_no_ : ww_pat
                        -> pattern
                     -> if &VAR_PGT^=_no_ : pgt
                        -> pattern
                     -> if &DATA_PST^=_no_: pst
                        -> pattern (nsimtest x)
                     -> if r>2: paircomp
          -> ausgabe -> if &DATA_PST^=_no_: ausgabe2
                     -> if r>2: ausgabe3
                     -> if &VAR_PIT^=_no_ : ausgabe4
                     -> if &VAR_PGT^=_no_ : ausgabe5

*/