%Lab5.1
%Clean,clean,clean!
close all
clear all

%GetValues
run L51values.m

%GetAPandTP'sPartialDerivatives
%RemoveFirstTwoPoints,PureH2OandEtOH
for int=start:last
    dens_dens_TP=dens_dens(int);
    temp_dens_TP=15.6;
    [x_Et,TP,AP]=trueproof(dens_dens_TP,temp_dens_TP);
    trial1(int)={[x_Et,TP,AP]};
end

%GetTrueProofandMolFracFromCells
%PutIntoMatrix
for int=start:last
    holder=[trial1{1,int}];
    mol_graph1(int)=holder(1);
    TP_graph1(int)=holder(2);
end

%GetAPandTP'sPartialDerivatives
%RemoveFirstTwoPoints,PureH2OandEtOH
for int=start:last
    dens_dens_TP=dens_dens(int);
    temp_dens_TP=20;
    [x_Et,TP,AP]=trueproof(dens_dens_TP,temp_dens_TP);
    trial2(int)={[x_Et,TP,AP]};
end

%GetTrueProofandMolFracFromCells
%PutIntoMatrix
for int=start:last
    holder=[trial2{1,int}];
    mol_graph2(int)=holder(1);
    TP_graph2(int)=holder(2);
end

%GetAPandTP'sPartialDerivatives
%RemoveFirstTwoPoints,PureH2OandEtOH
    for int=start:last
    dens_dens_TP=dens_dens(int);
    temp_dens_TP=25;
    [x_Et,TP,AP]=trueproof(dens_dens_TP,temp_dens_TP);
    trial3(int)={[x_Et,TP,AP]};
end

%GetTrueProofandMolFracFromCells
%PutIntoMatrix
for int=start:last
    holder=[trial3{1,int}];
    mol_graph3(int)=holder(1);
    TP_graph3(int)=holder(2);
end
