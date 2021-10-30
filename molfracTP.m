%MolFractoTrueProof
function[TP]=molfracTP(x_Et)    
    %GetValues
    run L51values.m
    
    %GetVolP_Et
    syms VolP_Et
    equ=(VolP_Et*Dens_Et/mm_Et)/...
        ((VolP_Et*Dens_Et/mm_Et)+((100-VolP_Et)*Dens_H2O/mm_H2O));
    VolP_solve=solve(equ==x_Et,VolP_Et);
      
    %GetTP
    TP=double(VolP_solve*2);
end