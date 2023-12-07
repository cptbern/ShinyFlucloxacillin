floxjointmap <- '
$PROB
# Joint model unbound (cmt=1) and total (cmt=2)

$PARAM @annotated
TVCL: 66.09 : Clearance typical
TVV: 217.39 : Volume typical
TVKB: 5.98 : binding rate albumin typical

ETA1: 0 : Volume IIV
ETA2: 0 : Clearance IIV
ETA3: 0 : KB IIV
ETA4: 0 : Volume IOV
ETA5: 0 : Clearance IOV
ETA6: 0 : KB IOV

$PARAM @annotated @covariates 
GFR: 70 :Kidney function
ALB: 20 : Albumin
WT: 70 : Weight
PROP_u: 0.0676 : proportional error total
PROP_t: 0.0144 : proportional error unbound
ADD_u: 0.0036 :additive error unbound
ADD_t: 0.7225 : additive error total

$CMT @annotated
CENT : Central() [ADM, OBS] // unbound
PERI: Peripheral() [OBS] // total

$SET end=100, delta=0.5

$MAIN
double VCi = TVV * exp(ETA(1)+ ETA1) * exp(ETA(4) + ETA4);
double CLi = TVCL*pow((GFR/70),1.09) *exp(ETA(2)+ ETA2) * exp(ETA(5) + ETA5);
double KBi = TVKB * pow((ALB/20),1.63) * exp(ETA(3)+ETA3) * exp(ETA(6) + ETA6);

$OMEGA @annotated
EV: 0.3025 : omega Volume IIV
ECL: 0.5476 : omega Clearance IIV
EKB: 0.2401 : omega KB IIV

GV: 0.2025 : gamma Volume IOV
GCL: 0.0961 : gamma Clearance IOV
GKB: 0.004096 : gamma KB IOV

$SIGMA
0 // proportional unbound
1 // additive unbound
0 // proportional bound
1 // additive bound


$ODE
dxdt_CENT = -(CLi/VCi)*CENT;

$TABLE
double IPRED_u = CENT/VCi;
double PAR = IPRED_u + (ADD_u + PROP_u*IPRED_u)*EPS(2)+EPS(1);

double IPRED = IPRED_u * KBi + IPRED_u;
double MET = IPRED + (ADD_t + PROP_t*IPRED)*EPS(4) + EPS(3);

double DV = PAR;
if(self.cmt == 2) DV = MET;

$CAPTURE DV PAR MET VCi CLi KBi
'
