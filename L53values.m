%Lab5.3:ValuesfromNotebook
%PowerValues
pow_boil=[0.61,0.53];
pow_chill=[2.94,3.38];
kwhtoJ=3.6*10^6;

%Density,Flow,MoleFractionsandTempValues
%ForProduct,Waste,andFeed
prod_flo=[5.6515,4.7772];
prod_dens=[0.8371,0.8252];
prod_temp=[23.3,22.4];
inp_dens=0.9808;
inp_temp=20.1;
x_wast=0.0201;
wast_flo=[128,41.225];
wast_dens=[0.9801,0.9767];
temp_in=[22.1,20.5,22.8,18.4,...
18.6,21.4,21.9,21.7];
dens_in=[0.81,0.8146,0.8142,0.8236,...
0.8301,0.8328,0.846,0.8727];
H2O_temp=[14.2,22.5];
cool_flo=[350,350];

%RefluxRatio
O=3;

%UnitConversions
ppm=1/10^6;
kilo=10^3;

%MolarMassandMassValues
mm_EtOH=46.07;
mm_H2O=18.06;
step=0.001;
L_mol=0:step:1;
start=1;
last=length(prod_dens);

%Antoine's Coefficients
%GammaCoefficients
A12=1.6789;
A21=0.9227;
A_H2O=7.96681;
B_H2O=1668.21;
C_H2O=228;
A_EtOH=8.04494;
B_EtOH=1554.30;
C_EtOH=222.650;
P_mmHg=760;

%LatentHeatCoefficients
Hvap_H2O=2258;
Hvap_EtOH=838;
last2=size(L_mol,2);