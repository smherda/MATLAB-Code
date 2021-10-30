%dpdproofanddproofdTFunction
function[dp_dprf,dprf_dT]=partialfrac(AP,Temp_pyc)
    %GetValues
    run L51values.m

    %GetTempinFarenheitGetNeighboringRows
    Temp_pyc_f=Temp_pyc*5/9+32;
    Temp_min=floor(Temp_pyc_f);
    Temp_max=ceil(Temp_pyc_f);

    %GetNeighboringApparentProofs
    AP_min=floor(AP);
    AP_max=ceil(AP);
        
    %GetDifferenceinProofFromUsingAPandTempValues
    dprf_dT=(TTB1(AP_max,Temp_max)-TTB1(AP_min,Temp_min))/...
        ((Temp_max-Temp_min)*9/5-64);
            
    %GetChangesBetweenTheseValuesAndDivideToGetAnswer
    dp_dprf=(TTB6(AP_max)-TTB6(AP_min))*Dens_H2O./(AP_max-AP_min)/...
        TTB1(AP_max,Temp_max)-TTB1(AP_min,Temp_min);
end