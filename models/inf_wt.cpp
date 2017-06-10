[PARAM] @annotated
CL        : 3.83    : Clearance (L/h)
V         : 39.4    : Volume (L) 
WT        : 70     : Weight (kg)

[FIXED] @annotated 
STD_WT   : 70  : Standard weight (kg)

[CMT] @annotated
CENT : Central compartment (mg)

[PKMODEL]
ncmt=1, trans=11

[MAIN]
D_CENT = 1; 

  
double CLi = CL*pow(WT/STD_WT,0.75)*exp(ECL);
double Vi = V*(WT/STD_WT)*exp(EV);

[OMEGA] @annotated @block @correlation @name IIV
ECL : 0.044       : Eta on CL
EV  : 0.66 0.04 : Eta on V

[SIGMA] @annotated
PROP : 0.04 : Proportional error (variance)
ADD  : 0.8  : Additive ((mg/L)^2)
  
[TABLE]
double CP = CENT/Vi;
double DV = CP*(1+PROP) + ADD;

[CAPTURE] @annotated
CP     : predicted plasma concentration (mg/L)
DV     : plasma concentration (mg/L)
CLi    : Individual Clearance (L/hr)
Vi     : Individual Volume (L)
WT     : Weight (kg)
