$PROB
Population PK model of imatinib in GIST patients, as published by Demetri et al.

$PARAM @annotated
TVCL :  8.18  : Clearance (L/hr)
TVV  :  168   : Volume of distribution (L)
KA   :  1     : Absorption rate (dummy)
F1   :  0     : Fraction zero-order absorption (dummy)
D1   :  1.69  : Duration of zero-order absorption (hr)

$PARAM @covariates @annotated
ALB  :  38    : Typical albumin value
ALBCL:  1.66  : Effect of albumin on CL
ALBV :  1.66  : Effect of albumin on V
WBC  :  7.76  : Typical WBC count
WBCCL: -0.418 : Effect of WBC on CL
WBCV : -0.418 : Effect of WBC on V

$OMEGA @correlation @block @annotated

ECL : 0.120                  : Random effect on CL
EV  : 0.119 0.128            : Random effect on V

$MAIN
D_CENT = D1;
F_CENT = F1;
F_GUT = 1-F1;

double CL = TVCL * pow((ALB/38.3), ALBCL) * pow((WBC/(7)), WBCCL) * exp(ECL);
double V = TVV * pow((ALB/38), ALBV) * pow((WBC/(7)), WBCV) * exp(EV);

$PKMODEL cmt = "GUT CENT", depot = TRUE                                                                            
                                                                            
$SIGMA @annotated

PROP: 0.122  : Proportional residual error
ADD : 1.58E-5   : Additive residual error

$TABLE
double IPRED = CENT / V;
double DV = IPRED + sqrt((pow(EPS(1),2)*(pow((CENT/V),2)) + (pow(EPS(2),2))));

$CAPTURE @annotated

IPRED   : Predicted plasma concentration (mg/L)
  DV      : Observed plasma concentration (mg/L)