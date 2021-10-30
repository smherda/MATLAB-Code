%Lab5.2:Calculations2
%Clean,clean,clean!
close all
clear all

%GetValues
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

%OrganizeforGraphs
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

%Power
power=[0.99,1.08,1.12,1.16,1.21,1.26];
p1=power(1);
p2=p1+power(2);
p3=p2+power(3);
p4=p3+power(4);
p5=p4+power(5);
power_tot=[p1,p2,p3,p4,p5];
power_tot_joul=power_tot.*3.6*10^6;

%InputLatentHeatandHeatCapacities
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

%InputConversions
P_mmHg=760;
approx=10^2-1;
ppm=0.01/10^6;
p=10;

%OtherSetValues
start=1;
last=5;

%GetInstantaneousMolFract
for int=start:last
    [x_vial_out]=trueproof(dens_vial(int),temp_vial(int)); 
    x_vial(int)=x_vial_out;
end

%GetFlaskMolFrac
for int=start:last
    [x_flask_out]=trueproof(dens_flask(int),temp_flask(int)); 
    x_flask(int)=x_flask_out;
end

%ConvertforEasierNotation
x1f=x_flask(1);
x2f=x_flask(2);
x3f=x_flask(3);
x4f=x_flask(4);
x5f=x_flask(5);
x1v=x_vial(1);
x2v=x_vial(2);
x3v=x_vial(3);
x4v=x_vial(4);
x5v=x_vial(5);

%GetVariablesforRKMethod
for int=1:approx
    %GetGammas
    g1=exp(A12*(A21*(1-x_EtOH)/...
        (A12*x_EtOH + A21*(1-x_EtOH))^2));
    g2=exp(A21*(A12*(x_EtOH)/...
        (A12*x_EtOH + A21*(1-x_EtOH)))^2);

    %UseGammatoGetTemperature  
    solv=@(temp_in)(x_EtOH*g1*10^(A_EtOH-B_EtOH/...
        (temp_in+C_EtOH))+(1-x_EtOH)*...
        g2*10^(A_H2O-B_H2O/(temp_in+C_H2O))-P_mmHg);
    
    %GetVaporMolFraction
    y_in=(x_EtOH*g1*10^(A_EtOH-B_EtOH/((fzero(solv,90))+C_EtOH)));
    y_in=y_in/P_mmHg;
    
    %GetX,Y,andLvalues
    h=-1/10;
    f=@(L_tot,x_in)((y_in-x_in)*(L_tot)^(-1));    
    d1=f(L_tot,x_EtOH);
    d2=f(L_tot+h/2,x_EtOH+d1*h/2);
    d3=f(L_tot+h/2,x_EtOH+d2*h/2);
    d4=f(L_tot+h,x_EtOH+h*d3);
    
    %CalculateEthanolandTotalMoles
    x_EtOH=x_EtOH+(d1+2*d2+2*d3+d4)*h/6;
    L_tot=h+L_tot;
       
    % Store the x, y, T, and L values
    xsave(int+1)=x_EtOH;
    ysave(int+1)=y_in;
end

%GetTheoreticalValues
x_pot_theo=xsave;  
x_flask_theo=ysave;

%GetPlotsforTheoreticalValues
potaxis=1:length(x_pot_theo);
flaskaxis=1:length(ysave);

%ErrorforErrorBarsfromManual
%GetUncertaintyinMolFrac
%FromtheManual
err_val=0.001;
err=ones(length(err_val))*err_val;
err=err(1,:);
dens_flask_err=err+dens_flask;

for int=start:last
    [x_flask_out]=trueproof(dens_flask_err(int),temp_flask(int)); 
    x_flask_err(int)=x_flask_out;
end

%ArrangeErrorsForGraph
err=ones(length(x_flask))*std(abs(x_flask-x_flask_err));
err=err(1,:);

%PlotTheoreticalandMeasuredValues
plot(flaskaxis,x_flask_theo)
hold on;
plot(potaxis,x_pot_theo)
hold on;

%LabelAxisandLegend
xlabel('Cumulative Mass (g)');
ylabel('Mole Fraction of Ethanol');
legend('Condensor, Theoretical','Pot, Theoretical');
hold on

%InputHistogramfortheFlasks
rectangle('Position',[0,0,m1f,x1f]);
rectangle('Position',[m1f,0,mass_flask(2),x2f]);
rectangle('Position',[m2f,0,mass_flask(3),x3f]);
rectangle('Position',[m3f,0,mass_flask(4),x4f]);
rectangle('Position',[m4f,0,mass_flask(5),x5f]);
alpha(0)
hold on

%InputHistogramfortheVials
rectangle('Position',[0,0,m1v,x1v]);
rectangle('Position',[m1v,0,mass_vial(2),x2v]);
rectangle('Position',[m2v,0,mass_vial(3),x3v]);
rectangle('Position',[m3v,0,mass_vial(4),x4v]);
rectangle('Position',[m4v,0,mass_vial(5),x5v]);
alpha(0)
hold on

%InputErrorbars
err_x=[m1f/2,m1f+mass_flask(2)/2,m2f+mass_flask(3)/2,...
    m3f+mass_flask(4)/2,m4f+mass_flask(5)/2];
x_flask=[x1f,x2f,x3f,x4f,x5f];
errorbar(err_x,x_flask,err,'LineStyle','none');
hold on

%GetValuesfromExcelSheet
area_vial=[1.33E+06,1.31E+06,...
    1.31E+06,1.01E+05,5.39E+05]*10^-3;
area_flask=[1.09E+06,2.18E+06,3.42E+06,...
    4.28E+06,2.42E+06]*10^-3;
dens_vial=[0.9715,0.977,0.9842,0.9926,0.9936];
dens_flask=[0.884,0.899,0.94,0.976,0.9936];

%GetvialantaneousMolFrac
for int=start:last
    [x_vial_out]=trueproof(dens_vial(int),temp_vial(int)); 
    x_vial(int)=x_vial_out;
end

%GetFlaskMolFrac
for int=start:last
    [x_flask_out]=trueproof(dens_flask(int),temp_flask(int)); 
    x_flask(int)=x_flask_out;
end

%LabelAxisandLegend
figure
xlabel('Mole Fraction of Ethanol');
ylabel('A/p (Millimeter per mole)');
hold on 
legend('Vials');
hold on
%InsertZeroforStart
x_plot=[x_vial(1,:),0];

%PlotFigure
y_vial=area_vial.*p.*(x_vial.*mm_EtOH+(1-x_vial).*mm_H2O)./dens_vial;
y_vial=[y_vial(1,:),0];
plot(x_plot,y_vial,'o');
hold on

%PlotHistogram
y_flask=area_flask.*(x_flask*mm_EtOH+(1-x_flask)*mm_H2O)./dens_flask;
rectangle('Position',[0,0,x_flask(5),y_flask(1)]);
rectangle('Position',[x_flask(5),0,x_flask(4)-x_flask(5),y_flask(2)]);
rectangle('Position',[x_flask(4),0,x_flask(3)-x_flask(4),y_flask(3)]);
rectangle('Position',[x_flask(3),0,x_flask(2)-x_flask(3),y_flask(4)]);
rectangle('Position',[x_flask(2),0,x_flask(1)-x_flask(2),y_flask(5)]);
alpha(0)
hold on

%GetTheoreticalforKValues
x_lin=linspace(0,x_vial(1));
y=polyval(polyfit(x_plot,y_vial,2),x_lin);

%GetRelationshipBetweenXandY
y2=polyval(polyfit(x_flask,y_flask,2),x_lin);
K=y./y2;
x_pot_theo=xsave*ppm;  
x_flask_theo=x_pot_theo.*K;

%PlotTheoreticalandMeasuredValues
figure
plot(flaskaxis,x_flask_theo)
hold on;
plot(potaxis,x_pot_theo)
hold on;

%LabelAxisandLegend
xlabel('Cumulative Mass (g)');
ylabel('Mole Fraction of O-Cresol');
legend('Condensor, Theoretical','Pot, Theoretical');
hold on
Ks=mean(K);

%InputHistogramfortheFlasks
rectangle('Position',[0,0,m1f,x1v*ppm*Ks]);
rectangle('Position',[m1f,0,mass_flask(2),x2v*ppm*Ks]);
rectangle('Position',[m2f,0,mass_flask(3),x3v*ppm*Ks]);
rectangle('Position',[m3f,0,mass_flask(4),x4v*ppm*Ks]);
rectangle('Position',[m4f,0,mass_flask(5),x5v*ppm*Ks]);
alpha(0)
hold on

%InputHistogramfortheVials
rectangle('Position',[0,0,m1v,x1v*ppm]);
rectangle('Position',[m1v,0,mass_vial(2),x2v*ppm]);
rectangle('Position',[m2v,0,mass_vial(3),x3v*ppm]);
rectangle('Position',[m3v,0,mass_vial(4),x4v*ppm]);
rectangle('Position',[m4v,0,mass_vial(5),x5v*ppm]);
alpha(0)
hold on

%Lab5.2:CalculationandPlot4
figure
plot(mass_flask_graph,power_tot_joul,'o');
hold on

%LabelAxisandLegend
xlabel('Mass Distilled (g)');
ylabel('Energy Consumed (J)');
legend('Energy, Measured','Energy, Theoretical');
hold on 

%UsetheHeatEquationtoGetthePower
temp_flask_EtOH=boil_EtOH-temp_flask;
temp_flask_H2O=boil_H2O-temp_flask;
pow_theo2=mass_flask.*x_flask*Cp_EtOH.*temp_flask_EtOH+...
    mass_flask.*(1-x_flask)*Cp_EtOH.*temp_flask_H2O;
pow_theo=mass_flask.*x_flask*Hvap_EtOH+mass_flask.*(1-x_flask)*Hvap_H2O;

%ArrangethePowerfortheGraph
%Cumulative,SoAddOntoEachOther
%TakeIntoAccountthe15PercentEfficiency
p1=pow_theo(1)+pow_theo2(1);
p2=pow_theo(2)+pow_theo2(2);
p3=pow_theo(3)+pow_theo2(3);
p4=pow_theo(4)+pow_theo2(4);
p5=pow_theo(5)+pow_theo2(5);
pow_theo_graph=[p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5]/0.15;
plot(mass_flask_graph,pow_theo_graph);

%GetLinearfitsforErrorAnalysis
polyfit(mass_flask_graph,pow_theo_graph,1);
polyfit(mass_flask_graph,power_tot_joul,1);