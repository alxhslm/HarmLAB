function X = time2freq(x,NHarm,Nfft)
X = fft(x,Nfft,1)/Nfft;
X = fftshift(X,1);

iz = floor(Nfft/2)+1;
X = X(iz + (0:NHarm),:);

X(2:end,:) = X(2:end,:) * 2;  %2 as one-sided