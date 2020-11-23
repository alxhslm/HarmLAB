function frf = hbm_excitation_forces(problem,frf)
for i = 1:length(frf)
    frf(i).Fe = (problem.Ku*frf(i).U.' + problem.Cu*(1i*frf(i).W.*frf(i).U).' + problem.Mu*(frf(i).W.^2.*frf(i).U).').';
end