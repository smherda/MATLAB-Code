%Lab5.2:ValuesfromNotebook1
%Temperatures
temp_flask=[22.3,21.1,23.1,24.8,24];
temp_vial=[24.2,26.3,28.3,27.4,31.3,30.2];
temp_col=[37.1,88.5,92.2,97,99,99.5];

%SetValuesofMass
mass_H2Os=344.195;
mass_EtOHs=67.793;
mm_H2O=18.01;
mm_EtOH=46.07;

%MolesandMolFrac
mol_H2Os=mass_H2Os/mm_H2O;
mol_EtOHs=mass_EtOHs/mm_EtOH;
L_tot=mol_EtOHs+mol_H2Os;
x_EtOH=mol_EtOHs/L_tot;

%MassValuesfromNotebookandCalculations
massV_full=[12.247,12.398,12.229,12.356,12.302,12.652];
massF_full=[132.192,135.966,134.996,135.09,130.35];
massV_emp=[6.67,6.685,6.695,6.688,6.669,6.683];
massF_emp=[95.507,96.859,94.581,94.183,92.522];
mass_flask=massF_full-massF_emp;
mass_vial=massV_full-massV_emp;

%OrganizeMassforUpcomingGraphs
m1f=mass_flask(1);
m2f=m1f+mass_flask(2);
m3f=m2f+mass_flask(3);
m4f=m3f+mass_flask(4);
m5f=m4f+mass_flask(5);
m1v=mass_vial(1);
m2v=m1v+mass_vial(2);
m3v=m2v+mass_vial(3);
m4v=m3v+mass_vial(4);
m5v=m4v+mass_vial(5);
mass_flask_graph=[m1f,m2f,m3f,m4f,m5f];

%Densities
dens_flask=[0.884,0.899,0.9842,0.9926,0.9936];
dens_vial=[0.9715,0.977,0.9842,0.9926,0.9936,0.9935];
dens_EtOH=0.79313;
dens_H2O=0.99904;

%VolCalculations
vol_pot=400;
vol_flask=mass_flask./dens_flask;
vol_vial=mass_vial./dens_vial;

%Power
power=[0.99,1.08,1.12,1.16,1.21,1.26];
p1=power(1);
p2=p1+power(2);
p3=p2+power(3);
p4=p3+power(4);
p5=p4+power(5);
power_tot=[p1,p2,p3,p4,p5];
power_tot_joul=power_tot.*3.6*10^6;

%LatentHeatandHeatCapacities
Hvap_H2O=2258;
Hvap_EtOH=838;
Cp_H2O=70;
Cp_EtOH=112;
boil_EtOH=79;
boil_H2O=100;

%Antoine'sCoefficientsandvanLaar'sParameters
A_H2O=7.96681;
B_H2O=1668.21;
C_H2O=228;
A_EtOH=8.04494;
B_EtOH=1554.30;
C_EtOH=222.650;
A12=1.6789;
A21=0.9227;

%Conversions
P_mmHg=760;
approx=10^2-1;
ppm=1/10^6;

%OtherSetValues
start=1;
last=5;