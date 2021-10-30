%ECH145A:Lab5,Tables
%InputSetValues
mm_Et=46.07;
mm_H2O=18.06;
Dens_Et=0.79313;
Dens_H2O=0.99904;
Cels_60=15.556;
alpha=25*10^-6;
start=1;
last=10;

%GetDatafromProvidedTables
TTB6=xlsread('TTB_Table_6_digitized.xlsx','E7:E208');
TTB1=xlsread('TTB_Table_1_digitized.xlsx','B8:CX208');

%GetDensityandTempsfromLab51
temp_dens=xlsread('Lab51.xlsx','L6:L15');
dens_dens=xlsread('Lab51.xlsx','K6:K15');
temp_pyc=xlsread('Lab51.xlsx','J6:J15');

%GetValuesforGraph
dens_pyc=dens_dens.*(1+alpha.*(temp_pyc-Cels_60));

%GetError
err_pyc_temp=0.5/mean(temp_pyc);
err_dens_dens=0.001;
error=sqrt(err_pyc_temp^2+err_dens_dens^2)/2;