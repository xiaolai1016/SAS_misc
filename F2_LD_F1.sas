/*********************************************************************
     Makro          F2_LD_F1

     Version 1.1 (DATUM Februar 2003)

     Die Variable VAR enth"alt die Werte
     Die Faktoren FACTOR1 und FACTOR2 enthalten die WHOLE-PLOT-F.
     Die Variable time enthält den SUB-PLOT-FACTOR
     der Faktor SUBJECT soll die Patienten/ Probanden enthalten
     QUELLE: - DFG, Br 655/11-1
             - Akritas, Brunner Langer ....
     geschrieben von Steffen Ballerstedt und
                     Andreas Oelerich, AMS Goettingen, FEB  1998
     letzte Version : 04.03.2003 Carola Werner
*********************************************************************/

%macro F2_LD_F1(DATA= ,VAR = ,FACTOR1 = ,
                FACTOR2 = ,TIME = ,SUBJECT = );

/********************************************************************/
%MACRO middle(dat=,var=,fak1=,fak2=,pat=,t=);

 PROC SORT DATA=&dat;
  BY &t &pat &fak2 &fak1;

 PROC RANK DATA = &dat TIES = mean out = _daten_;
  VAR &var;
  RANKS rijk;

 PROC FREQ DATA = _daten_;
  TABLES &fak1*&fak2*&pat*&t / OUT = _ni_ NOPRINT;

 PROC FREQ DATA = _ni_;
  TABLES &fak1          / OUT=ni1 NOPRINT;

 PROC FREQ DATA = _ni_;
  TABLES &fak2          / OUT=ni2 NOPRINT;

 PROC FREQ DATA = _ni_;
  TABLES &t             / OUT=ni3 NOPRINT;


 PROC FREQ DATA = _ni_;
  TABLES &fak1*&fak2    / OUT=ni12 NOPRINT;

 PROC FREQ DATA = _ni_;
  TABLES &fak1*&t       / OUT=ni13 NOPRINT;

 PROC FREQ DATA = _ni_;
  TABLES &fak2*&t       / OUT=ni23 NOPRINT;

 PROC FREQ DATA = _ni_;
  TABLES &fak1*&fak2*&t / OUT=ni123 NOPRINT;

RUN;

 PROC SORT DATA = _daten_;
  BY &fak1 &fak2 &pat &t;

 PROC MEANS DATA = _daten_ NOPRINT;
  BY &fak1 &fak2 &pat &t;
  VAR rijk;
  OUTPUT OUT = _value_ N = mijk MEAN = RmIJK_;

 DATA _value_;
  SET _value_;
  KEEP mIJK RmIJK_ &fak1 &fak2 &pat &t;

 RUN;
%MEND;

%middle(dat=&data,var=&var,fak1=&factor1, fak2=&factor2, pat=&subject,t=&time);

proc iml worksize=120;

RESET LINESIZE=80;

/**********************************************************************/
START cov (NN,RIJK_,ni,d,V);
/***********************************************************************
    Die Funktion cov berechnet die Kovarianzmatrix des

                          general mixed models

    Die vorliegende Funktion berechnet die Kovarianz des multivariaten
    Modells, Spezialfälle wie das CS-Modell müssen extra programmiert
    werden.

    QUELLE: Vorlesungsskript
    'Theory of rank Tests in factorial designs' von Prof. E. Brunner
    15. Oktober 1997.
    Programmiert von:       Andreas Oelerich, AMS Göttingen, Januar 1998

    Die Parameter sind
    NN    : Anzahl der observation
    RIJK_ : die Ränge wobei über die Messwiederholung gemittelt ist
            (MACRO middle!)
    ni    : Der Vektor der die Anzahl der subjects enthält
    d     : Anzahl der abhängigen Levels (repeated measures)
    V     : Die Kovarianz die zurückgegeben wird
***********************************************************************/
n = SUM(ni);                        /* Anzahl der subjects            */
r = NROW(ni);                       /* Anzahl der unabhängigen Levels */
pijk_ = 1/NN * (RIJK_ - 1/2);       /* Die im Vektor RIJK_ enthaltenen
                                       Ränge sollen die Mittelränge
                                       über die Messwiederholungen der
                                       subjects sein -> MACRO middle!
                                       pijk_ sind die Phi_ijk. des
                                       Papers                         */
IND = LOC(pijk_=.);                 /* wo sind missings?              */
IF NCOL(IND) >0 THEN pijk_[IND] = 0;
pijk_ = SHAPE(pijk_,n,d);           /* ... als Matrix                 */
lijk  = J(n,d,1);                   /* die lambda_ijk des Papers      */
IF NCOL(IND) > 0 THEN lijk[IND] = 0;/* Kodierung !                    */
lp    = lijk # pijk_ ;
DO i=1 TO r;
   IF i=1 THEN DO; r1 = 1;                  r2 = ni[i];        END;
          ELSE DO; r1 = SUM(ni[1:i-1]) + 1; r2 = SUM(ni[1:i]); END;
   lij_   = (lijk[r1:r2,])[+,];     /* ist Vektor!!                   */
   l_1lp  = lp[r1:r2,];
   pij__  = (1 / lij_) # (l_1lp[+,]);
                                    /* Phi_ij..                       */
   vij_d  = J(d,d,0);
   vijj_d = J(d,d,0);
   DO j=1 TO d;
      pij_h = J(ni[i],1,1) @ pij__[j]; /* pij_h ist eine Hilfsmatrix
                                       die die Mittel (pij__)
                                       aufplustert, so ist die
                                       Operation "-" unten zulässig
                                       wegen gleicher Dimension       */
      vij_d[j,j] = SUM(( (lijk[r1:r2,j]) # (pijk_[r1:r2,j] - pij_h))##2) ;
      faktor2    = ni[i] / (lij_[j] * ( lij_[j] - 1 ) );
      vij_d[j,j] = vij_d[j,j] * faktor2;
                                    /* vij_d sind die Schätzer der
                                       Diagonalelemente der Kovarianz-
                                       matrix Vn                      */
      DO j2=1 TO d;
         IF j2 ^= j
            THEN DO;
                   pij_h2 = J(ni[i],1,1) @ pij__[j2];
                   lijj_        = SUM(lijk[r1:r2,j] # lijk[r1:r2,j2]);
                   vijj_d[j,j2] = SUM((lijk[r1:r2,j]) # (lijk[r1:r2,j2])
                    # (pijk_[r1:r2,j] - pij_h) # (pijk_[r1:r2,j2] - pij_h2) );
                   faktor       = ni[i] / ((lij_[j]-1) * (lij_[j2]-1)
                                  + lijj_ - 1);
                   vijj_d[j,j2] = faktor * vijj_d[j,j2];
                 END;               /* of THEN DO                     */
      END;                          /* of DO : j2                     */
   END;                             /* of DO : j                      */
   Vi = n / ni[i] * (vij_d + vijj_d);
   IF i=1 THEN V = Vi;
          ELSE V = BLOCK(V,Vi);
END;                                /* of DO : i                      */
print V;
FINISH;                             /* of FUNCTION : cov              */
/**********************************************************************/
start wald(NN,C,V,p,df,Q,p_val);
cvc= C * V * C`;
g_cvc = ginv(cvc);
Q     = NN * (C * p)` * (g_cvc) * (C * p);
df=trace(cvc * g_cvc);
if Q >= 0 then p_val = 1 - probchi(Q,df);
else p_val = .;
finish;

START box(N,p_dach,C,V_d,B1,fg,p);
T        = C` * GINV(C * C`) * C;
tv       = T * V_d;
B1       = (N * p_dach` * T * p_dach) / (trace(tv));
fg       = ((trace(tv))**2)/(trace(tv * tv));
fg2      = 10000;
if B1 >= 0 & fg > 0 then p = 1 - probf(B1,fg,fg2);
else p = .;
FINISH;
/********************************************************************/
/* box2 berechnet nach der Reduzierung (Zeit rausmitteln) lediglich */
/* die geschaetzten Freiheitsgrade (gemaess fixed model) von den    */
/* 3 Tests  fuer  A , B , AB                                        */
/********************************************************************/
start box2(r,a,b,c,ni,N,v_d,fg_box2);

 lambda = inv(diag(ni)-I(r));
 mat    = I(a) @ I(b) @ j(1,c,1/c);

 /* s_d ist reduzierte Kovarianzmatrix mit Diagonalgestalt */
 s_d    = mat * v_d * mat`;
 pa = (I(a)-j(a,a,1/a));  da = diag(pa @ j(b,b,1/b));
 pb = (I(b)-j(b,b,1/b));  db = diag(j(a,a,1/a) @ pb);
                         dab = diag(pa @ pb);


 fg_box2=j(3,1,0);
 fg_box2[1]=((trace(da *s_d))**2)/(trace(da *da *s_d*s_d*lambda));
 fg_box2[2]=((trace(db *s_d))**2)/(trace(db *db *s_d*s_d*lambda));
 fg_box2[3]=((trace(dab*s_d))**2)/(trace(dab*dab*s_d*s_d*lambda));

 *B1 = (N*p_dach`*T*p_dach)/trace(tv);
 *fg = ((trace(tv))**2)/(trace(tv*tv));

finish;  /* box2 */


/*********************************************************************
   Einlesen der Rangmittel und des Designs
*********************************************************************/
USE _value_;
READ ALL VAR{mijk}     INTO m_ijk;
READ ALL VAR{RmIJK_}   INTO RmIJK;
READ ALL VAR{&factor1} INTO fak1_;
READ ALL VAR{&factor2} INTO fak2_;
READ ALL VAR{&subject} INTO pat_;
READ ALL VAR{&time}    INTO t_;
CLOSE _value_;

ranks = RmIJK;
lev_a = unique(fak1_);                /* Die Stufen des Faktors A ***/
lev_b = unique(fak2_);                /* Die Stufen des Faktors B ***/
lev_c = UNIQUE(t_);                   /* Die Stufen des Faktors C ***/
lev_d = unique(pat_);                 /* Die Stufen des Faktors D ***/
a     = ncol(lev_a);                  /* Anzahl der Stufen von A  ***/
b     = ncol(lev_b);                  /* Anzahl der Stufen von B  ***/
c     = ncol(lev_c);                  /* Anzahl der Stufen von C  ***/
ni = J(a * b,1,0);

DO i=1 TO a;
   DO j=1 TO b;
      IND    = LOC((fak1_=lev_a[i]) # (fak2_ = lev_b[j]));
      ni[(i-1)*b+j] = NCOL(UNIQUE(pat_[IND]));
   END;
END;

nsum = SUM(ni);

IF nsum*c ^= nrow(m_ijk)
  THEN DO;
       print 'Error',
             'If you have missing values (empty cell(s))',
             'please fill in one observation with a .';
       END;
m_ijk = SHAPE(m_ijk,nsum,c);
N     = sum(m_ijk);                   /* Anzahl der Daten           */

/********************************************************************/
/*** i ist der Laufindex des Faktors A  i=1,...,na (whole1)       ***/
/*** j ist der Laufindex des Faktors B  j=1,...,nb (whole2)       ***/
/*** k ist der Laufindex des Faktors C  k=1,...,nc (sub)          ***/
/*** l ist der Laufindex der Subjects   l=1,...,nd                ***/
/*** s ist der Laufindex der Messwiederholungen s=1,...,m_ijk     ***/
/********************************************************************/
 p_a    = I(a)-j(a,a,1/a) ;
 p_b    = I(b)-j(b,b,1/b) ;
 p_c    = I(c)-j(c,c,1/c) ;

 eins_a =   j(1,a,1/a) ;
 eins_b =   j(1,b,1/b) ;
 eins_c =   j(1,c,1/c) ;

 ca    = p_a     @  eins_b  @  eins_c ;
 cb    = eins_a  @  p_b     @  eins_c ;
 cc    = eins_a  @  eins_b  @  p_c    ;
 cab   = p_a     @  p_b     @  eins_c ;
 cac   = p_a     @  eins_b  @  p_c    ;
 cbc   = eins_a  @  p_b     @  p_c    ;
 cabc  = p_a     @  p_b     @  p_c    ;

 ta    = ca`  * ginv(ca  *ca` ) * ca ;
 tb    = cb`  * ginv(cb  *cb` ) * cb ;
 tc    = cc`  * ginv(cc  *cc` ) * cc ;
 tab   = cab` * ginv(cab *cab`) * cab;
 tac   = cac` * ginv(cac *cac`) * cac;
 tbc   = cbc` * ginv(cbc *cbc`) * cbc;
 tabc  = cabc`* ginv(cabc*cabc`)* cabc;

/********************************************************************/
/*** Rijk    = Vector of rank means for cell (i,j,k)              ***/
/*** Rijk_ma = (anxbc)-Matrix der Rangmittelwerte                 ***/
/*** ranks   = R_ijks                                             ***/
/*** p_d     = estimate for relativ treatment effect              ***/
/********************************************************************/
p       = (ranks - 1/2) / N ;
r_ijk   = ranks;       /* Matrix der Rangmittel ueber die Mess-     */
                       /* wiederholungen                            */
r_ijk   = shape(r_ijk,nsum,c);
/********************************************************************/
/*** In der Matrix r_ijk stehen die Rangmittel der Patienten      ***/
/********************************************************************/
r_ij  = j(a*b,c,0);
m_ij  = j(a*b,c,0);
do i=1 to a*b;
   if i>1 then do; r1 = sum(ni[1:i-1]) + 1; end;
          else do; r1 = 1;                  end;
   r2 = sum(ni[1:i]);
   r_ij[i,] = (r_ijk[r1:r2,])[:,];
   m_ij[i,] = (m_ijk[r1:r2,])[+,];
end;
p_s = 1/N * (r_ij - 1/2);
p_d = shape(p_s,a*b*c,1);

RUN COV (N,RmIJK,ni,c,V);

V = V * N / nsum;

/*********************************************************************

     Die WALD-TYP Statistk

*********************************************************************/

Q     = j(7,1,0);
p_val = j(7,1,0);
df    = j(7,1,0);
/*df[1] = a-1;
df[2] = b-1;
df[3] = c-1;
df[4] = (a-1) * (b-1);
df[5] = (a-1) * (c-1);
df[6] = (b-1) * (c-1);
df[7] = (a-1) * (b-1) * (c-1);*/

/*********************************************************************
     Wald-Typ Statistiken
*********************************************************************/

run wald(N,CA   ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt A     */
df[1]=dfi; Q[1]  = QF; p_val[1]  = p_v;
run wald(N,CB   ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt B     */
df[2]=dfi; Q[2]  = QF; p_val[2]  = p_v;
run wald(N,CC   ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt C     */
df[3]=dfi; Q[3]  = QF; p_val[3]  = p_v;
run wald(N,CAB  ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt AB    */
df[4]=dfi; Q[4]  = QF; p_val[4]  = p_v;
run wald(N,CAC  ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt AC    */
df[5]=dfi; Q[5]  = QF; p_val[5]  = p_v;
run wald(N,CBC  ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt BC    */
df[6]=dfi; Q[6]  = QF; p_val[6]  = p_v;
run wald(N,CABC ,V ,p_d,dfi ,QF,p_v) ;            /* Effekt ABC   */
df[7]=dfi; Q[7]  = QF; p_val[7]  = p_v;

wald = Q || df || p_val;

/*********************************************************************
   Die BOX-Approximation
*********************************************************************/
Q_box = j(7,1,0);
p_box = j(7,1,0);
dfbox = j(7,1,0);
RUN box(N,p_d,CA   ,V,QF_,DF_,p_);             /* Effekt A         */
Q_box[1]  = QF_; dfbox[1]  = DF_ ; p_box[1]  = p_;
RUN box(N,p_d,CB   ,V,QF_,DF_,p_);             /* Effekt B         */
Q_box[2]  = QF_; dfbox[2]  = DF_ ; p_box[2]  = p_;
RUN box(N,p_d,CC   ,V,QF_,DF_,p_);             /* Effekt C         */
Q_box[3]  = QF_; dfbox[3]  = DF_ ; p_box[3]  = p_;
RUN box(N,p_d,CAB  ,V,QF_,DF_,p_);             /* Effekt AB        */
Q_box[4]  = QF_; dfbox[4]  = DF_ ; p_box[4]  = p_;
RUN box(N,p_d,CAC  ,V,QF_,DF_,p_);             /* Effekt AC        */
Q_box[5]  = QF_; dfbox[5]  = DF_ ; p_box[5]  = p_;
RUN box(N,p_d,CBC  ,V,QF_,DF_,p_);             /* Effekt BC        */
Q_box[6]  = QF_; dfbox[6]  = DF_ ; p_box[6]  = p_;
RUN box(N,p_d,CABC ,V,QF_,DF_,p_);             /* Effekt ABC       */
Q_box[7]  = QF_; dfbox[7]  = DF_ ; p_box[7]  = p_;

ANOVA = Q_box || Dfbox || p_box ;

/*********************************************************************
   Die BOX2 - Approximation, schätzen des 2. Freiheitsgrades
*********************************************************************/
RUN box2(a*b,a,b,c,ni,N,V,dfbox2);
box2 = j(3,4,0);  /* 3Zeilen fuer A,B,AB, 4Spalten fuer T,fg1,fg2,p */
box2[,3] = dfbox2;             /* fg2 fuer F(fg1,fg2) - fuer A,B,AB */
box2[1,2]= dfbox[1];           /* fg1 fuer F(fg1,fg2) - fuer A      */
box2[2,2]= dfbox[2];           /* fg1 fuer F(fg1,fg2) - fuer B      */
box2[3,2]= dfbox[4];           /* fg1 fuer F(fg1,fg2) - fuer AB     */
box2[1,1]= Q_box[1];           /* Statistik fuer A                  */
box2[2,1]= Q_box[2];           /* Statistik fuer B                  */
box2[3,1]= Q_box[4];           /* Statistik fuer AB                 */
if box2[1,1] >= 0 & box2[1,2] > 0 &  box2[1,3] > 0
  then box2[1,4] = 1 - probf(box2[1,1], box2[1,2], box2[1,3]);
else box2[1,4] = .;
if box2[2,1] >= 0 & box2[2,2] > 0 &  box2[2,3] > 0
  then box2[2,4] = 1 - probf(box2[2,1], box2[2,2], box2[2,3]);
else box2[2,4] = .;
if box2[3,1] >= 0 & box2[3,2] > 0 &  box2[3,3] > 0
  then box2[3,4] = 1 - probf(box2[3,1], box2[3,2], box2[3,3]);
else box2[3,4] = .;
/********************************************************************/
/*** Diese Funktion schreibt die Faktorstufenkombinationen zweier   */
/*** Faktoren in den Vektor help.                                   */
/********************************************************************/
start out1(ch_1,ch_2,ch1_n,ch2_n,help);
do i=1 to ch1_n;
   do j=1 to ch2_n;
      help = help // (ch_1[i] + '*' + ch_2[j]);
   end;
end;
finish; /* out1                                                     */
/********************************************************************/
/*** Diese Funktion schreibt die Faktorstufenkombinationen dreier   */
/*** Faktoren in den Vektor help                                    */
start out2(ch_1,ch_2,ch_3,ch1_n,ch2_n,ch3_n,help);
 do i=1 to ch1_n;
    do j=1 to ch2_n;
       do k=1 to ch3_n;
          help = help // (ch_1[i] + '*' + ch_2[j] + '*' + ch_3[k]);
       end;
    end;
 end;
finish; /* out2                                                     */
/********************************************************************/
/*** Diese Funktion erstellt den source-Vektor der fuer den Output  */
/*** benoetigt wird.                                                */
/********************************************************************/
start source_k(a,b,c,cha_1,cha_2,cha_3,source) ;
%let name1 = %trim(&factor1);
%let name2 = %trim(&factor2);
%let name3 = %trim(&time);
spalte1 =j(a,1,"&name1") // j(b,1,"&name2") // j(c,1,"&name3")
                         // j(a*b,1,"&name1"+'*'+"&name2"+' ')
                         // j(a*c,1,"&name1"+'*'+"&name3"+' ')
                         // j(b*c,1,"&name2"+'*'+"&name3"+' ')
                         // j(a*b*c,1,"&name1"+'*'+"&name2"+'*'+"&name3"+' ');
run out1(cha_1,cha_2,a,b,text);
run out1(cha_1,cha_3,a,c,text);
run out1(cha_2,cha_3,b,c,text);
run out2(cha_1,cha_2,cha_3,a,b,c,text);
spalte2 = (cha_1`) // (cha_2`) // (cha_3`) // text;
source  = concat (spalte1,spalte2);
finish; /* source_k                                                 */
/********************************************************************/
/********************************************************************/
/*** Ausgabe der Quadratformen und der p_values                  ****/
/********************************************************************/
start test_out(wald,anova,box2);
source =  {"A" , "B" , "T" , "AB" , "AT" , "BT" , "ABT" ,
           "A|BT" , "B|AT" , "T|AB"};
row    =  {"A" , "B" , "T" , "AB" , "AT" , "BT" , "ABT"};
row2   =  {"A", "B" , "AB"};
col1   = {"QF"  "DF"  "p-value"};
col2   = {"QF", "DF1" , "DF2" , "p-value"};
wald   = wald[1:7,];
anova  = anova[1:7,];
print /;
print 'Wald-type statistic',
      'Approximation for large sample sizes with Chi^2_DF';
print wald[format=6.5][r=row][c=col1];
print 'ANOVA-type statistic',
      'Box-approximation for small sample size with Chi^2_(df)',
      '-estimated degree of freedom';
print anova[format=6.5][r=row][c=col1];
print 'Anova-type-statistic',
      'modified Box-Approximation for the whole-plot factors A and B'
      'for small sample size with F(df1,df2)',
      '-estimated degrees of freedom';
print box2[format=6.5][r=row2][c=col2];
finish;
/********************************************************************/
start CLI(a,b,c,nsum,m_ijk);
 reset center;
 class  = {A &factor1, B &factor2, T &time};
 levels = a // b // c ;
 print / ;
 print 'F2_LD_F1 --- subjects(A x B) x T',
       'A, B, T: fixed, subjects: random';
 print "SAS-datafile-name: &data",
       "Response variable: &var";
 print 'Class Level Information' ;
 print class levels;
 reset noname;
 print 'Total number of observations ' (sum(m_ijk)) ' ',
       'Total number of subjects     ' nsum ' ',
       'Number of missing values     ' (NCOL(LOC(m_ijk=0)))  ;
finish ; /* Ende von O_CLI */;
/********************************************************************/
/********************************************************************/
start nu_char(v,v_neu) ;
if type (v) ='N' then v_neu = char(v,max(int (log10(max(abs(v))))+1));
 else v_neu = trim(v);
finish;
/* Ende der Funktion nu_char                                       **/
/********************************************************************/

start rmeans(a,b,c,source,mijk,p_d,r_ijk);
m_ijk  = shape(mijk,a*b*c,1);
N = SUM(mijk);
p_ijk = p_d;
r_ijk = shape(r_ijk,a*b*c,1);
/*CB    = j(1,a*c,1) @ i(b);
CC    = j(1,a,1) @ i(c) @ j(1,b,1) ;
CAB   = i(a) @ j(1,c,1) @ i(b);
CAC   = i(a) @ i(c) @ j(1,b,1);
USE ni1  ; READ ALL VAR{count} INTO n_a  ; CLOSE ni1  ;
USE ni2  ; READ ALL VAR{count} INTO n_b  ; CLOSE ni2  ;
USE ni3  ; READ ALL VAR{count} INTO n_c  ; CLOSE ni3  ;
USE ni12 ; READ ALL VAR{count} INTO n_ab ; CLOSE ni12 ;
USE ni13 ; READ ALL VAR{count} INTO n_ac ; CLOSE ni13 ;
USE ni23 ; READ ALL VAR{count} INTO n_bc ; CLOSE ni23 ;
USE ni123; READ ALL VAR{count} INTO n_abc; CLOSE ni123;*/

CA    = I(a) @ j(1,b*c,1);
CB    = J(1,a,1) @ I(b) @ J(1,c,1);
CC    = J(1,a*b,1) @ I(c);
CAB   = I(a*b) @ J(1,c,1);
CAC   = I(a) @ j(1,b,1) @ I(c);
CBC   = j(1,a,1) @ I(b*c);

/*nobs*/
n_a   = CA  * m_ijk;
n_b   = CB  * m_ijk;
n_c   = CC  * m_ijk;
n_ab  = CAB * m_ijk;
n_ac  = CAC * m_ijk;
n_bc  = CBC * m_ijk;
n_abc = m_ijk;

nor = n_a // n_b // n_c // n_ab // n_ac // n_bc // n_abc;

/*rankmeans*/
r_ijkx=r_ijk#m_ijk;
r_a   = CA  * r_ijkx / n_a;
r_b   = CB  * r_ijkx / n_b;
r_c   = CC  * r_ijkx / n_c;
r_ab  = CAB * r_ijkx / n_ab;
r_ac  = CAC * r_ijkx / n_ac;
r_bc  = CBC * r_ijkx / n_bc;
r_abc = r_ijk;


r_vec = r_a // r_b // r_c // r_ab // r_ac // r_bc // r_abc ;
p = (r_vec -0.5)/N;
print 'RTE  = Relative Treatment Effects          ',
      'Nobs = Number of observations (do not count',
      'the repeated measurements within the cells)';
print source[c={source}] r_vec[format=6.5 c={"Rank mean"}]
      nor[c={"Nobs"}] p[format=6.5 c={"RTE"}];
finish;  /* Ende von rang_vec */
/*******************************************************************/
/*** Ausgabe der Warnungen                                       ***/
/*******************************************************************/
start warning(c,n,V,m_ijk);
if c > n | det(V) < 10**(-10) then
  print / '';
reset nocenter;
if c > n then
  print 'Warning:',
        'There are less subjects than sub-plot factor levels.';
if det(V)< 10**(-10) then
  print 'Warning:',
        'Do not use the Wald-type-statistic, because the covariance matrix is singular.';
*if min(eigval(V)) < 0 & NCOL(LOC(m_ijk=0)) > 0 then
*  print 'Warning:',
*        'The estimated covariance matrix is not positive semidefinite due to missing values.';
reset center;
finish;


/********************************************************************/
dataname = name(&data);
varname  = name(&var);
fa1_name = name(&factor1);
fa2_name = name(&factor2);
fa3_name = name(&time);
fa4_name = name(&subject);
run nu_char(lev_a, leva);
run nu_char(lev_b, levb);
run nu_char(lev_c, levc);
run nu_char(lev_d, levd);
run cli(a,b,c,nsum,m_ijk);
run source_k(a,b,c,leva,levb,levc,source);
RUN rmeans(a,b,c,source,m_ij,p_d,r_ij);
RUN warning(c,n,V,m_ijk);
RUN test_out(wald,anova,box2);

/********************************************************************/
quit; /* IML */
/********************************************************************/
%mend F2_LD_F1;
/********************************************************************/




