%Lab5.1:Plots
%GetData
run L51calc.m

%GetValuesforNewDensity
for int=start:last
    [x_Et,TP,AP]=trueproof(dens_pyc(int),temp_pyc(int));
    [p1,p2]=partialfrac(AP,temp_pyc(int));
    p1store(int)=p1;
    p2store(int)=p2;
end

%ConvertValuesforDensity
dens_new=dens_pyc'+p1store.*p2store.*(temp_pyc-temp_dens)';

%RemoveN/APoints
dens_dens(2)=[];
dens_pyc(2)=[];
dens_new(2)=[];
mol_graph1(2)=[];
mol_graph2(2)=[];
mol_graph3(2)=[];

%PlotandLabel,fig1
xlabel('Density, Pycnometer (g/cm^3)');
ylabel('Density, Densitometer (g/cm^3)');
xlim([min(dens_pyc),max(dens_pyc)]);
hold on
plot(dens_pyc,dens_dens,'o');
hold on
plot(dens_pyc,dens_dens);
hold on
herrorbar(dens_pyc,dens_dens,error*ones(1,length(dens_dens)));
hold on
errorbar(dens_pyc,dens_dens,err_dens_dens*ones(1,length(dens_dens)));
hold on

%GetEquation
polysave=polyfit(dens_pyc,dens_dens,1);
text(0.9,0.86,['y = ',num2str(polysave(1)),'+',num2str(polysave(2))]);

%InputFiguresfortheSecondGraph
figure
hold on

%InputTheoreticalLines
plot(sort(dens_dens,'descend'),sort(mol_graph1));
hold on
plot(sort(dens_dens,'descend'),sort(mol_graph2));
hold on
plot(sort(dens_dens,'descend'),sort(mol_graph3));
hold on

%InputPoints
plot(dens_dens,mol_graph1,'o');
hold on
plot(dens_dens,mol_graph2,'o');
hold on
plot(dens_dens,mol_graph3,'o');
hold on

%LabeltheLinesontheGraph
xlabel('Density, Densitometer (g/cm^3)');
ylabel('Mol Fraction');
yyaxis right
ylabel('True Proof');
legend('15.6C','20C','25C');

%SetUpTicksforRightAxis
tick=0:0.1:1;
for int=1:11
    TP_graph(int)=molfracTP(tick(int));
end
yticks(tick);

%LabeltheTicks
yticklabels({num2str(TP_graph(int-10)),num2str(TP_graph(int-9)),...
num2str(TP_graph(int-8)),num2str(TP_graph(int-7)),...
num2str(TP_graph(int-6)),num2str(TP_graph(int-5)),...
num2str(TP_graph(int-4)),num2str(TP_graph(int-3)),...
num2str(TP_graph(int-2)),num2str(TP_graph(int-1)),...
num2str(TP_graph(int))});
