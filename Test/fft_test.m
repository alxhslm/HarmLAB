function fft_test 

NHarm = [4 4]; 
Nfft = [128 128];
iSub = [2 7;
        4 1
        2 4
        3 1];
for i = 1:10000
    X = rand(4,1) + rand(4,1)*1i;
    tic;
    x = freq2time3d(X,NHarm,iSub,Nfft,1);
    t(i) = toc;
    tic;
    x = freq2time3d(X,NHarm,iSub,Nfft,0);
    t2(i) = toc;
end
mean(t)
mean(t2)

function x = freq2time3d(X,NHarm,iSub,Nfft,b)
%NB, scaling so IFFT = Expanding Fourier Coefficients

%first regenerate the -ve frequency components
Xfft = unpackfreq(X,NHarm,iSub);

if b
    %now augment for Nfft frequencies
    iz = floor(Nfft/2)+1;
    Xfft2 = zeros(Nfft(1),Nfft(2),size(X,3));
    Xfft2(iz(1)+(-NHarm(1):NHarm(1)),iz(2)+(-NHarm(2):NHarm(2)),:) = Xfft;
    Xfft = Xfft2;
end

%finally compute the ifft
Xfft = ifftshift(Xfft,1);
Xfft = ifft(Xfft,Nfft(1),1)*Nfft(1);
if NHarm(2)>0 
    Xfft = ifftshift(Xfft,2);
    x = ifft(Xfft,Nfft(2),2)*Nfft(2);
else
    x = Xfft;
end

%the answer should now be real to machine precision
x = real(x);

%wrap hypertime
x = reshape(x,prod(Nfft),[]);



function x = freq2time(X,NHarm,Nfft,b)

%first regenerate the -ve frequency components
X = cat(1,0*X(2:end,:,:), X);
Xfft = 0.5*(X + conj(flip(X,1)));

if b
    %now augment for Nfft frequencies
    iz = floor(Nfft/2)+1;
    sz = size(X);
    Xfft2  = zeros([Nfft,sz(2:end)]);
    Xfft2(iz+(-NHarm:NHarm),:,:) = Xfft;
    Xfft = Xfft2;
end

%finally compute the ifft
Xfft = ifftshift(Xfft,1);
x    = ifft(Xfft,Nfft,1) * Nfft;

%the answer should now be real to machine precision
x    = real(x);