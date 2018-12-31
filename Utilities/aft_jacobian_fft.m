function Ju = aft_jacobian_fft(harm)
Nfft = harm.Nfft;
NComp = harm.NComp;

theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
[theta1,theta2] = ndgrid(theta1,theta2);
theta0 = [theta1(:)'; theta2(:)'];

theta1 = permute(theta0(1,:),[1 3 2]);
theta2 = permute(theta0(2,:),[1 3 2]);

ju = zeros(NComp,1,Nfft(1)*Nfft(2));
ju(1,1,:) = 1/(Nfft(1)*Nfft(2));
for l = 1:(harm.NFreq-1)
    ph = harm.kHarm(l+1,1)*theta1 + harm.kHarm(l+1,2)*theta2;
    cl = cos(ph); sl = sin(ph);
    ju(2*l,:,:)   =  2*cl/(Nfft(1)*Nfft(2));
    ju(2*l+1,:,:) = -2*sl/(Nfft(1)*Nfft(2));
end
Ju = ju;