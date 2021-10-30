%Lab5:TrueProofFunction
function[x_Et,TP,AP]=trueproof(Dens,Temp)
     %GetValues
     run L51values.m

    %TakeDiffToFindClosestRows
    SG_dens=Dens/Dens_H2O;
    SG_diff=abs(SG_dens-TTB6);
    C2r_low=find(min(SG_diff)==SG_diff);
    C2r_high=C2r_low+1;

    %CalculateApparentProof
    diff=(Dens-TTB6(C2r_low))./...
        (TTB6(C2r_high)-TTB6(C2r_low))*(C2r_high-C2r_low);
    AP=C2r_low+diff;   
    AP_high=ceil(AP);
    AP_low=floor(AP);
    
    %ConvertCelsiustoFarenheit 
    TempF_low=floor(Temp*9/5+32);
    TempF_high=ceil(Temp*9/5+32);
    
    %InterpolateforTrueProof
    %GetPartialFractions
    dfdC=TTB1(AP_high,TempF_low)-TTB1(AP_low,TempF_low);
    dfdT=TTB1(AP_low,TempF_high)-TTB1(AP_low,TempF_low);
    TP1=TTB1(AP_low,TempF_low)+(AP-AP_high)*(dfdC);
    TP2=(TempF_high-TempF_low)*(dfdT);
    TP=TP1+TP2;
    
    %GetVolumePercent
    volper_Et=TP*0.5;   
    volper_H2O=abs(volper_Et-100);
    
    %GetMolFracofEthanol
    per_Et=volper_Et*Dens_Et;
    num=per_Et/mm_Et;
    div=per_Et/mm_Et+(volper_H2O*Dens_H2O/mm_H2O);
    x_Et=num/div;
end