%Lab5.2:Calculations1
%Clean,clean,clean!
close all
clear all

%GetValues
run L52values1.m

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

%GuessforTempandFzero
guess=90;

%GetVariablesforRKMethod
for int=start:approx
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
    y_in=(x_EtOH*g1*10^(A_EtOH-B_EtOH/((fzero(solv,guess))+C_EtOH)));
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
potaxis=start:length(x_pot_theo);
flaskaxis=start:length(ysave);

%ErrorforErrorBarsfromManual
%GetUncertaintyinMolFrac
%FromtheManual
err_val=0.001;
err=ones(length(err_val))*err_val;
err=err(1,:);
dens_flask_err=err+dens_flask;

%FindErroratEachPoint
for int=start:last
    D_flask_in=dens_flask_err(int);
    T_flask_in=temp_flask(int);
    [x_flask_out]=trueproof(D_flask_in,T_flask_in); 
    x_flask_err(int)=x_flask_out;
end

%PutErrorIntoStandardFormationforGraphing
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
