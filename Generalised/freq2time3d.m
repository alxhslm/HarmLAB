function x = freq2time3d(X,NHarm,iSub,Nfft)
%NB, scaling so IFFT = Expanding Fourier Coefficients

%first regenerate the -ve frequency components
X = unpackfreq(X,NHarm,iSub);

%now augment for Nfft frequencies
iz = floor(Nfft/2)+1;
Xfft = zeros(Nfft(1),Nfft(2),size(X,3));

for i = 1:2
    iHarm{i} = (-NHarm(i):NHarm(i));
    if NHarm(i) == 0
        iHarm{i}  = 1;
    end
end
Xfft(iz(1)+iHarm{1},iz(2)+iHarm{2},:) = X;

%finally compute the ifft
Xfft = ifftshift(Xfft,1);
Xfft = ifft(Xfft,Nfft(1),1)*Nfft(1);

Xfft = ifftshift(Xfft,2);
x = ifft(Xfft,Nfft(2),2)*Nfft(2);


%the answer should now be real to machine precision
x = real(x);

%wrap hypertime
x = reshape(x,prod(Nfft),[]);
