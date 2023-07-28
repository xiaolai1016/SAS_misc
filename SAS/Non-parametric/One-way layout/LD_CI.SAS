/*********************************************************************************************
* Filename: LD_CI.SAS
* Author:   Sebastian Domhof
* Source:   Brunner, Domhof, Langer: Nonparametric Analysis of Longitudinal Data
* Date:     26. Oct. 2000
* letzte Version : 04.03.2003 Carola Werner
**********************************************************************************************
* Description:
* This SAS-Macro computes
* - estimators for the relative treatment effects
* - estimators for the bias of the previous estimators
* - confidence intervals for the relative treatment effects
* in a design with at most one whole-plot- and at most one sub-plot factor. The numbers of 
* subjects in each group may be different, but it is assumed that there is exactly one 
* observation for each subject at each level of the sub-plot-factor.
**********************************************************************************************
* Usage:
* %LD_CI(DATA=data, VAR=var, GROUP=group, TIME=time, SUBJECT=subject, ALPHA=alpha)
* with
* - data:    name of the SAS-data-set
*            default: _last_
* - var:     name of the variable containing the observations
* - group:   name of the variable containing the group(=whole-plot-factor-level)-names
*            if there is just one group, this variable may be left unspecified
* - time:    name of the variable containing the time(=sub-plot-factor-level)-names
*            if there is just one time, this variable may be left unspecified
* - subject: name of the variable containing the subject-names
* - alpha:   1 - the confidence level
*            default: 0.05
* All subjects must have different names. If every group starts with subject 1 again, this
* macro will crash!
*********************************************************************************************/

%MACRO LD_CI(DATA=_last_, VAR=, GROUP=_none_, TIME=_none_, SUBJECT=, ALPHA=0.05);

PROC IML;
RESET NONAME; *print no variable names;

/*********************************************************************************************
* Description of the variables: (in brackets: variable names in the source)
* - variables read directly from the SAS-data-set
*   var      (cv): observations (X_{iks})
*   grp_org  (cv): groups like in the original data-set
*   subj_org (cv): subjects like in the original data-set
*   time_org (cv): times like in the original data-set
*
* - variables for counting the observations
*   nn           : total number of observations (N)
*   n            : number of subjects (n)
*   a            : number of groups (a)
*   t            : number of time-points (t)
*   n_i      (cv): number of subjects per group (n_i)
*
* - variables for organizing groups, subjects, times
*   grp_nm   (cv): set of sorted group-names
*   group    (cv): groups with names replaced by numbers 1, ..., a
*   grp_out  (cv): group-numbers according to the vector of relative effects
*   subj_nm  (cv): set of subject-names
*   subject  (cv): groups with names replaced by numbers 1, ..., n_1, ..., 1, ..., n_a
*   time_nm  (cv): set of sorted time-point-names
*   time     (cv): groups with names replaced by numbers 1, ..., t
*   time_out (cv): time-point-numbers according to the vector of relative effects
*
*   matrix    (m): matrix to compute R_{ik.} from R_{iks} and R_{uk.}^{(-is)} from
*                  R_{ukw}^{(-is)}
*
* - variables containing different rankings of the observations
*   r_iks    (cv): overall ranks (R_{iks})
*   r_ik     (cv): sums of ranks over time (R_{ik.})
*   r_row_i  (cv): used for indexing r_ik with group-numbers
*   r_iks_i  (cv): group-internal ranks (R_{iks}^{(i)})
*   r_iks_ik (cv): subject-internal ranks (R_{iks}^{(ik)})
*   r_iks_is (cv): ranks within the group-time-combinations (R_{iks}^{(is)})
*   r_ukw_is  (m): ranks without group-time-combinations (R_{ukw}^{(-is)}), ukw=row,
*                  is=column, (R_{iks}^{(-is)} = 0 by definition)
*   r_uk_is   (m): sums of ranks without group-time-combinations over time (R_{uk.}^{(-is)}),
*                  uk=row, is=column
*
* - estimators:
*   p_is     (cv): relative effects (p_{is})
*   b_s_i    (cv): bias (b_s^{(i)})
*   p_is_s   (cv): linearly transformed relative effects (p_{is}^*)
*   sigma_is (cv): variances (\sigma_{is}^2)
*   cl        (m): confidence-limits (p_{is,L}, p_{is,U}
*********************************************************************************************/

/*********************************************************************************************
* Sort and read the data-set. If there is just one group or time, generate the group or time
* variable, respectively. 
*********************************************************************************************/

IF "&TIME"^="_none_" THEN
  SORT DATA=&DATA BY &TIME;
SORT DATA=&DATA BY &SUBJECT;
IF "&GROUP"^="_none_" THEN
  SORT DATA=&DATA BY &GROUP;

USE &DATA;
READ ALL VAR{&SUBJECT} INTO subj_org;
READ ALL VAR{&VAR} INTO var;
IF "&GROUP"^="_none_" THEN
  READ ALL VAR{&GROUP} INTO grp_org;
ELSE
  grp_org=J(NROW(var),1,'_none_');
IF "&TIME"^="_none_" THEN
  READ ALL VAR{&TIME} INTO time_org;
ELSE
  time_org=J(NROW(var),1,'_none_');
CLOSE &DATA;


/*********************************************************************************************
* Initialize the values and dimensions of the variables
*********************************************************************************************/

* counting variables;
*********************************************************************************************;
nn=NROW(var);               * number of observations;

n=NCOL(UNION(subj_org));    * number of subjects;
a=NCOL(UNION(grp_org));     * number of groups;
t=NCOL(UNION(time_org));    * number of time-points;
n_i=J(a,1);                 * number of subjects per group;


* organization of groups, subjects, times;
*********************************************************************************************;
grp_nm=UNION(grp_org)`;     * set of group-names;
group=J(nn,1);              * group-numbers in var, r_iks, ...;
grp_out=J(a*t,1);           * group-numbers in p_is, b_s_i, ...;
subj_nm=UNION(subj_org)`;   * set of subject-names;
subject=J(nn,1);            * subject numbers in var (restarting with 1 for each group);
time_nm=UNION(time_org)`;   * set of time-point-names;
time=J(nn,1);               * time-numbers in var, r_iks, ...;
time_out=J(a*t,1);          * time-point-numbers in p_is, b_s_i, ...;

matrix=J(n,n*t);            * matrix with r_ik=matrix*r_iks, r_uk_is=matrix*r_ukw_is;

* rankings;
*********************************************************************************************;
r_row_i=J(n,1);             * group-numbers in r_ik;
r_iks_i=J(nn,1);            * group-internal ranks;
r_iks_ik=J(nn,1);           * subject-internal ranks;
r_iks_is=J(nn,1);           * ranks within group-time-combinations;

* estimators;
*********************************************************************************************;
p_is=J(a*t,1);              * relative effects;
b_s_i=J(a*t,1);             * bias;
p_is_s=J(a*t,1);            * transformed relative effects;
sigma_is=J(a*t,1);          * variances;
cl=J(a*t,2);                * confidence-limits;


/*********************************************************************************************
* Computations concerning the layout of the data-set
*********************************************************************************************/

* check if there is exactly one observations per subject and time-point;
*********************************************************************************************;
DO k=1 TO n;
  DO s=1 TO t;
    IF NCOL(LOC(subj_org=subj_nm[k] & time_org=time_nm[s]))=0 THEN DO;
      PRINT "No observation at subject" (subj_nm[k]) "and time-point" (time_nm[s]) "!";
      ABORT;
    END;
    ELSE IF NCOL(LOC(subj_org=subj_nm[k] & time_org=time_nm[s]))>1 THEN DO;
      PRINT "More than one observation at subject" (subj_nm[k])
            "and time-point" (time_nm[s]) "!";
      ABORT;
    END;
	IF var[(k-1)*t+s]=. then do; PRINT "The macro does not support missing values. ";
      ABORT; 
    END;
  END;
END;

* initialize group, time, subject;
*********************************************************************************************;
DO iks=1 TO nn;
  group[iks]=LOC(grp_nm=grp_org[iks]);
  time[iks]=LOC(time_nm=time_org[iks]);
  IF iks=1 THEN                               * first subjects gets number 1;
    subject[iks]=1;
  ELSE IF group[iks]^=group[iks-1] THEN       * every group starts with subject 1;
    subject[iks]=1;
  ELSE IF subj_org[iks]^=subj_org[iks-1] THEN * next subject gets next number;
    subject[iks]=subject[iks-1]+1;
  ELSE                                        * same subject gets same number;
    subject[iks]=subject[iks-1];
END;

* initialize n_i, grp_out, time_out;
*********************************************************************************************;
row=1;
DO i=1 TO a;
  n_i[i]=MAX(subject[LOC(group=i)]);
  Do k=1 TO n_i[i];
    r_row_i[row]=i;
    row=row+1;
  END;
  DO s=1 TO t;
    grp_out[(i-1)*t+s]=i;
    time_out[(i-1)*t+s]=s;
  END;
END;
grp_out=grp_nm[grp_out];
time_out=time_nm[time_out];

* initialize matrix;
*********************************************************************************************;
matrix=I(n)@J(1,t);


/*********************************************************************************************
* Compute the rankings
*********************************************************************************************/

* initialize r_iks, r_ik;
*********************************************************************************************;
r_iks=RANKTIE(var);
r_ik=matrix*r_iks;

* initialize r_iks_i, r_iks_ik, r_iks_is, r_ukw_is, r_uk_is;
*********************************************************************************************;
r_ukw_is=REPEAT(var,1,a*t); max=MAX(var); row=1;
DO i=1 TO a;
  loc_i=LOC(group=i);
  r_iks_i[loc_i]=RANKTIE(var[loc_i]);
  DO k=1 TO n_i[i];
    loc_ik=LOC(group=i & subject=k);
    r_iks_ik[loc_ik]=RANKTIE(var[loc_ik]);
  END;
  DO s=1 TO t;
    loc_is=LOC(group=i & time=s);
    r_iks_is[loc_is]=RANKTIE(var[loc_is]);
    r_ukw_is[loc_is,(i-1)*t+s]=max+1; * to compute the ranks without observations from group;
                                      * i and time s, replace these obs. by max+1;
    r_ukw_is[,(i-1)*t+s]=RANKTIE(r_ukw_is[,(i-1)*t+s]);
    r_ukw_is[loc_is,(i-1)*t+s]=0;
  END;
END;
r_uk_is=matrix*r_ukw_is;


/*********************************************************************************************
* Compute the estimators
*********************************************************************************************/
DO i=1 TO a;
  loc_i=LOC(r_row_i=i);
  DO s=1 TO t;
    is=(i-1)*t+s;
    loc_is=LOC(group=i & time=s);

    * compute p_is;
    *****************************************************************************************;
    p_is[is]=(r_iks[loc_is][:]-1/2)/nn;

    * compute b_s_i;
    *****************************************************************************************;
    b_s_i[is]=n_i[i]*(r_iks_ik[loc_is][:]-1/2-(r_iks_i[loc_is][:]-1/2)/n_i[i])/nn/(n_i[i]-1);
	*print r_iks_ik[loc_is][:];

    * compute sigma_is;
    *****************************************************************************************;
    r_s_vec=2*r_iks[loc_is]-r_iks_is[loc_is]-r_ik[loc_i]+r_uk_is[loc_i,is];
    lambd_is=SSQ(r_s_vec-r_s_vec[:])/(n_i[i]-1)/nn**2;
    tau_uis=J(a,1,0);
    DO u=1 TO a;
      IF u^=i THEN DO;
        loc_u=LOC(r_row_i=u);
        r_s_vec=r_ik[loc_u]-r_uk_is[loc_u,is];
        tau_uis[u]=SSQ(r_s_vec-r_s_vec[:])/(n_i[u]-1)/n_i[i]**2;
      END;
    END;
    sigma_is[is]=n*lambd_is/n_i[i]+(n_i#tau_uis)[+]/n/t**2;

    * compute the confidence limits;
    *****************************************************************************************;
    p_is_s[is]=(nn*p_is[is]-n_i[i]/2)/(nn-n_i[i]);
    IF 0<p_is_s[is] & p_is_s[is]<1 THEN DO;
      p_isl_g=LOG(p_is_s[is]/(1-p_is_s[is]))
              -nn*SQRT(sigma_is[is]/n)*PROBIT(1-&ALPHA/2)/(nn-n_i[i])/p_is_s[is]
               /(1-p_is_s[is]);
      p_isu_g=LOG(p_is_s[is]/(1-p_is_s[is]))
              +nn*SQRT(sigma_is[is]/n)*PROBIT(1-&ALPHA/2)/(nn-n_i[i])/p_is_s[is]
               /(1-p_is_s[is]);
      cl[is,1]=n_i[i]/2/nn+(nn-n_i[i])*EXP(p_isl_g)/nn/(1+EXP(p_isl_g));
      cl[is,2]=n_i[i]/2/nn+(nn-n_i[i])*EXP(p_isu_g)/nn/(1+EXP(p_isu_g));
    END;
    ELSE DO;
      cl[(i-1)*t+s,1]=.;
      cl[(i-1)*t+s,2]=.;
    END;
  END;
END;

/*********************************************************************************************
* Print the output
*********************************************************************************************/

* data-set information;
*********************************************************************************************;
PRINT "LD_CI",
      "Bias-Estimation and Confidence-Intervals for Relative Effects";

PRINT "SAS-Data-Filename: &DATA",
      ({"Response-Variable:", "Group-Variable:", "Time-Variable:", "Subject-Variable"})
      ({"&VAR", "&GROUP", "&TIME", "&SUBJECT"})
      ({"     ", "     ", "     ", "     "})
      ({"Observations:", "Groups:", "Timepoints:", "Subjects:"})
      (nn // a // t // n);

* estimators for relative effects, biases and variances;
*********************************************************************************************;
PRINT "Relative Effects, Biases, Variances and Confidence-Limits (alpha=&ALPHA)";
colname1="Group";
colname2="Time";
colname3={"RE" "Bias" "Variance" "lower" "upper"};
IF "&GROUP"="_none_" THEN
  PRINT time_out[COLNAME=colname2]
        (p_is || b_s_i || sigma_is || cl)[COLNAME=colname3 FORMAT=6.5];
ELSE IF "&TIME"="_none_" THEN
  PRINT grp_out[COLNAME=colname1]
        (p_is || b_s_i || sigma_is || cl)[COLNAME=colname3 FORMAT=6.5];
ELSE
  PRINT grp_out[COLNAME=colname1] time_out[COLNAME=colname2]
        (p_is || b_s_i || sigma_is || cl)[COLNAME=colname3 FORMAT=6.5];

QUIT;

%MEND;
