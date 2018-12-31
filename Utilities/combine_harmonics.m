function Xfft2 = combine_harmonics(Xfft,harm)
NDof = sum([harm.group{:}.NDof]);
Xfft2 = zeros(harm.NFreq,NDof);
for i = 1:length(harm.group)
    Xfft2(harm.group{i}.iFreq,harm.group{i}.iDof) = Xfft{i};
end