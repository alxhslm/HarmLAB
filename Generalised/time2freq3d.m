function X = time2freq3d(x,NHarm,iSub,Nfft)

%unwrap the hypertime
x = reshape(x,Nfft(1),Nfft(2),[]);

%perform the actual fft operation
Xfft = fft(x,Nfft(1),1)/Nfft(1);
Xfft = fftshift(Xfft,1);

Xfft = fft(Xfft,Nfft(2),2)/Nfft(2);
Xfft = fftshift(Xfft,2);

%truncate to frequencies we want
iz = floor(Nfft/2)+1;
Xfft = Xfft(iz(1)+(-NHarm(1):NHarm(1)),iz(2)+(-NHarm(2):NHarm(2)),:);

%pack up the frequencies
X = packfreq(Xfft,iSub);
