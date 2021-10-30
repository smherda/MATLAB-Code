function [HL,HV]=enth(x,temp,y)
%VirialHeatCapacityCoefficientsforH2O
A_H2O=92.05;
B_H2O=3.995*10^-2;
C_H2O=-2.110*10^-4; 
D_H2O=5.347*10^-7;

%VirialHeatCapacityCoefficientsforEtOH
A_EtOH=59.34;
B_EtOH=3.635*10^-1;
C_EtOH=-1.216*10^-3;
D_EtOH=1.803*10^-6;
T0=298;

%HeatofFusion
Hf_H2O=39.50; 
Hf_EtOH=39.40; 
Hf_H2Os=6.0097; 
Hf_EtOHs=-277.38; 

%SpecificHeatsforH2OandEthanol
Cp_H2O=A_H2O*temp+(B_H2O/2)*temp^2+(C_H2O/3)*temp^3+(D_H2O/4)*temp^4;
Cp_EtOH=A_EtOH*temp+(B_EtOH/2)*temp^2+(C_EtOH/3)*temp^3+(D_EtOH/4)*temp^4;
Cp_H2Os=A_H2O*T0+(B_H2O/2)*T0^2+(C_H2O/3)*T0^3+(D_H2O/4)*T0^4;
Cp_EtOHs=A_EtOH*T0+(B_EtOH/2)*T0^2+(C_EtOH/3)*T0^3+(D_EtOH/4)*T0^4;
    
%GetEnthalpyofLiquidandVapor
HL=min((x)*(Cp_EtOH-Cp_EtOHs)+(1-x)*(Cp_H2O-Cp_H2Os));
HV=min(y*(Hf_EtOH-Hf_EtOHs+(Cp_EtOH-Cp_EtOHs))+...
    (1-y)*(Hf_H2O-Hf_H2Os+(Cp_H2O-Cp_H2Os))); 
end

