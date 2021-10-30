%Lab5.3:CalculationsandPlot
%Clean,clean,clean!
close all
clear

%GetValues
run L53values.m

for int3=start:last2
    %GetGammas
    x_find=(int3-1)*step;
    gam1=exp(A12*(A21*(1-x_find)/...
        (A12*x_find + A21*(1-x_find)))^2);
    gam2=exp(A21*(A12*(x_find)/...
        (A12*x_find + A21*(1-x_find)))^2);

    %UseGammatoGetTemperature  
    equ=@(T_in)(x_find*gam1*10^(A_EtOH-B_EtOH/...
        (T_in+C_EtOH))+(1-x_find)*gam2*...
        10^(A_H2O-B_H2O/(T_in+C_H2O))-P_mmHg);
    T_in=fzero(equ,90);

    %GetVaporMolFraction
    Psat_Et=10^(A_EtOH-B_EtOH/(T_in+C_EtOH));      
    y_out(int3)=(x_find*gam1*Psat_Et)/P_mmHg;
end

%PlotUpperandLowerLimitsforLine
plot(L_mol,y_out);
hold on
plot(L_mol,L_mol);
hold on

%ConvertMolesforOpperatingLines
for int=start:length(dens_in)
    [x_EtOH(int)]=trueproof(dens_in(int),temp_in(int));
end

%GetSlopeAndGraphLiquidOppLine
F=prod_flo(1);
y2=(O+F)/(O+1);
b2=x_wast*(1-F)/(O+1);

%GetIntoaLineforGraph
for int=start:length(dens_in)
    liq(int)=y2.*x_EtOH(int)+b2;
end

%PlotGraph
syms yz
polysav2=polyfit(x_EtOH,liq,1);
plo2=polysav2(1)*yz+polysav2(2);
fplot(plo2,'--');

%GetInformationforUpwardFlux
x_prod=trueproof(prod_dens(start),prod_temp(start));
y_p=O/(O+1);
b=x_prod/(O+1);

%CreateLoopForNewValuesforGraph
for int=start:length(dens_in)
    Upp_s(int)=x_EtOH(int)*y_p+b;
end

%PlotDataforUpwardFlux
syms yv
polysav=polyfit(x_EtOH,Upp_s,1);
plo=polysav(1)*yv+polysav(2);
fplot(plo,'--')
hold on

%LabelGraph
ylabel('Vapor Fraction, Ethanol');
xlabel('Liquid Fraction, Ethanol');
legend('Vapor-Liquid Equilibrium Line',...    
    'Equal Vapor and Liquid Mole Fractions Line',...
    'Lower Operating Line',...
    'Upper Operating Line');

%LimitsSinceMoleFractionNotAbove1
xlim([0,1]);
ylim([0,1]);

%DrewLinesonMatlabAfterwardsSince
%HavingTroubleInputtingStepsto
%RepresenttheNumberofPlates:[
%UsedToGetaGeneralIdeaOfTheNumberofPlates
annotation('line',[0.299846625766871,0.169478527607362],...
[0.379160636758321,0.377713458755427]);
annotation('line',[0.299079754601227,0.299079754601227],...
[0.555716353111433,0.379160636758321]);
annotation('line',[0.523773006134969,0.299079754601227],...
[0.554716353111431,0.554269175108536]);
annotation('line',[0.522239263803681,0.523006134969325],...
[0.648335745296671,0.554269175108538]);
annotation('line',[0.640337423312883,0.522239263803681],...
[0.646888567293777,0.646888567293777]);
annotation('line',[0.170245398773006 0.170245398773006],...
[0.379160636758321 0.183791606367583]);
annotation('line',[0.137269938650307 0.17101226993865],...
[0.185238784370478 0.185238784370476]);


%GetBoilerTrayOuput
for int=start:length(temp_in)
    [H_L(int),H_V(int)]=enth(x_EtOH(int),temp_in(int),x_EtOH(int));
end

%GettheHeatUsedBytheCondensor
Q_con=(O+1)*(H_V(length(temp_in))-H_V(1));

%FindCoefficienceofPerformanceandEfficiency
eff_chill=-(Q_con*kilo)/(pow_chill(1)*kwhtoJ);
COP=-(Q_con*kilo)/(pow_chill(1)*kwhtoJ+pow_boil(1)*kwhtoJ);

%O-CresolAmount
%TookMoleFractionsofEtOHforInput,Waste,andProd
%andFoundCorrespondingKinLab5.2
%GenerallyAroundFollowingNumber
K=3.5721;
x_oc_in=x_EtOH(start)*K*ppm;
x_oc_prod=x_prod*K*ppm;
x_oc_waste=x_EtOH(length(x_EtOH))*K*ppm;