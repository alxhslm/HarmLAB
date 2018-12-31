function U = test_excite(hbm,problem,w0)
P = problem.P;
U = zeros(hbm.harm.NFreq,1);
U(1) = P.f0;
U(2) = P.f;
if hbm.harm.NHarm(2)>0
    ii = hbm.harm.kHarm(:,1) == 0 & hbm.harm.kHarm(:,2) == 1; 
    U(ii) = P.f2;
end