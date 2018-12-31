function x = freq2time(X,NHarm,Nfft)
%NB, scaling so IFFT = Expanding Fourier Coefficients

%first regenerate the -ve frequency components
X = cat(1,0*X(2:end,:,:), X);
X = 0.5*(X + conj(flip(X,1)));

%now augment for Nfft frequencies
iz = floor(Nfft/2)+1;
sz = size(X);
Xfft  = zeros([Nfft,sz(2:end)]);
Xfft(iz+(-NHarm:NHarm),:,:) = X;

%finally compute the ifft
Xfft = ifftshift(Xfft,1);
x    = ifft(Xfft,Nfft,1) * Nfft;

%the answer should now be real to machine precision
x    = real(x);