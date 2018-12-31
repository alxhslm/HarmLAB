function Fe = hbm_excitation(hbm,problem,w0,U)
 w = hbm.harm.kHarm*(repmat(hbm.harm.rFreqBase',1,size(w0,2)) .* w0);
w = permute(w,[1 3 2]);
Wu = repmat(1i*w,1,size(U,2),1);
Udot  = U.*Wu;
Uddot = Udot.*Wu;

Fe = mtransposex(mtimesx(problem.Ku,mtransposex(U)) + mtimesx(problem.Cu,mtransposex(Udot)) + mtimesx(problem.Mu,mtransposex(Uddot)));