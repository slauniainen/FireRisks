
function [P,z,ETleafcrit]=tree_water_potential(ETmm, daylength,P0, L, LAI, k0s, h0s, a, bTs, bs)

%based on Couvreur et al Plant Cell and Environmnet for the simplified case of
%- linear change of P50 along the stem
%- Q, r and alpha contant along stem; this is obtained by setting Eleaf, Huber and Kx_sat as independent of z
%i.e., using eq 17 and 18
%in general, pressures appear as being reported with sign

%inputs
%ETmm: daily total transpiration per unit ground area (mm)
%daylength: hrs between sunrise and sunset (hrs)
%P0: water botential at the bottom of the canopy (MPa); can be assumed to be equal to soil water potential if neglecting soil to root conductance
%LAI: leaf area index (m2 leaf/m2 ground)
%L: hydraulic path length, approximated by canopy height (m)
%k0s: sdependece of saturated hydraulic conductivity at the base of the tree (z = 0) on L
%h0s: coefficients of the dependence of Huber on L; give Huber in cm2/m2
%a: parameter of the logistic vulnerability function; specifically a/4 is the slope of the vulnerability function at P50 (i.e., b) (1/MPa)
%bTs: parameter to determine bT, and hence b (or P50) i.e., the pressure at conductance equal to 1/Kx_sat, logistic vulnerability function; specifically, bT is the P50 at z=L (MPa); 
%bs: parameter to determine b (or P50) i.e., the pressure at conductance equal to 1/Kx_sat, logistic vulnerability function; specifically, bs is the slope of the b to z relationship

%internal parameters
dz=0.1; %spatial resolution (m)
rhow=10^3; %kg/m3
g=9.81; %m/s2
alpha=1; %dimensionless

%position along the hydraulix path (m)
z=[0:dz:L];

%pressure at 50% loss of conductance along the stem (MPa)
bT=bTs(1)+bTs(2)*L;
b=bT+bs*(L-z);

%slope of the dB/dz relation, assuming linear dependence of b on z
B=-bs+rhow*g/10^6*alpha; %MPa/m

%huber, i.e. Asap/Aleaf (assumed constant with z)
h0=h0s(1)+h0s(2)*L; %cm2sap/m2leaf
huber=h0/10^4; %10^4 to get m2 sap/m2leaf

%saturated hydraulic conductivity a
k0=k0s(1)+k0s(2)*L; %saturated hydraulic conductivity at the base of the tree (z = 0) (kg/(m s MPa))
Kx_sat=k0; %assuming kexp=0, i.e., constant k along the plant - NB: if setting kexp=0, there should be a divided by 2, but apparently this was somewhat incorporated in K0

%transformation of daily ET into subdaily values, with a focus on the time of the day at which the transpiration is highest (and the conditions likely more prone to fires)
%assume a parabolic dependence (see Mathematica file for the calculations)
ET=3*ETmm/10^3/(2*daylength*3600);  %transipiration rate, per unit ground area (m/s, i.e., m3water/m2ground/s)


%transform transpiration per unit ground area (m/s) in stem flux rate (kg/s)
ETleaf=ET/LAI; %m3H20/m2leaf/s
display('ETleaf in mmol/m2/s')
ETleafmmol=ETleaf*10^9/18 

%ETcrit, i.e. max achievable transpiration (in mmol/m2/s)
kexp=0;hs=0; %due to the focus on the simplified case
Z50=0.93; %immaterial for kexp=0
% ETthresh=1;
% if ETleafmmol>ETthresh %for very low Et, there is no need to determine ETcrit (and the code may become very slow)
%     ETcrit=ecritfnc(P0,L,a,bT,bs,k0,Z50,kexp,h0,hs);
% else
%     ETcrit=ETthresh;
% end
% 
%turn off the check on ETcrit
 ETcrit=100;

%check that E does not exceed ETcrit, where the model breaks donw
if ETleafmmol>ETcrit
    display('Transpiration exceeds max possible value')
    display('P(z) is determined assuming 95% of ETcrit, which corresponds to a fraction of actual ET')
    0.95*ETcrit/ETleafmmol
    ETleafcrit=0.95*ETcrit/10^9*18;
    ETleaf=ETleafcrit;
else
    ETleafcrit=NaN;
end


%resulting term Qr, to include in the formula for Y (Qr is in MPa/m)
Qr=ETleaf*rhow/Kx_sat/huber;

%inverse of the functional loss of conductivity, at each height z
b0=bT+bs*L; %P50 at the bottom of the canopy
Y0=exp(a*(P0-b0));
Y=Y0.*exp(-a*(Qr+B).*z)+B./(Qr+B)*(1-exp(-a*(Qr+B).*z));

%water potential, for each height z
P=b+ 1/a*log(Y-1);

return