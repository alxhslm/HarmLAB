function Xfft2 = split_harmonics(Xfft,harm)
for i = 1:length(harm.group)
    Xfft2{i} = Xfft(harm.group{i}.iFreq,harm.group{i}.iDof);
end