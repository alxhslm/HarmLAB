function O = hbm_output(hbm,problem,w0,u,x)
NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
Nfft  = hbm.harm.Nfft(1);
kHarm = hbm.harm.kHarm(:,1);
rBase = hbm.harm.rFreqBase(1);

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;

w = kHarm*rBase*w0;
t = (0:(Nfft-1))/Nfft*2*pi/(rBase*w0);

%now compute the fourier coefficients for xdot, xddot
X = unpackdof(x,NFreq-1,problem.NDof,iRetain);
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%convert to time domain
x      = freq2time(X,NFreq-1,Nfft);
xdot   = freq2time(Xdot,NFreq-1,Nfft);
xddot  = freq2time(Xddot,NFreq-1,Nfft);

U = unpackdof(u,NFreq-1,problem.NInput);
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

u     = freq2time(U,NFreq-1,Nfft);
udot  = freq2time(Udot,NFreq-1,Nfft);
uddot = freq2time(Uddot,NFreq-1,Nfft);

%push through the nl system
o = feval(problem.model,'output',t, x.',xdot.',xddot.',u.',udot.',uddot.',hbm,problem,w0).';

%finally convert into a fourier series
O = time2freq(o,NFreq-1,Nfft);

O = packdof(O);