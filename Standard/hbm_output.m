function O = hbm_output(hbm,problem,w0,u,x)
NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
Nfft  = hbm.harm.Nfft(1);
kHarm = hbm.harm.kHarm(:,1);
wBase = hbm.harm.rFreqBase(1)*w0(1);

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;

w = kHarm*wBase;
States.t = (0:(Nfft-1))/Nfft*2*pi/wBase;

%now compute the fourier coefficients for xdot, xddot
X = unpackdof(x,NFreq-1,problem.NDof,iRetain);
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%convert to time domain
States.x      = freq2time(X,NFreq-1,Nfft).';
States.xdot   = freq2time(Xdot,NFreq-1,Nfft).';
States.xddot  = freq2time(Xddot,NFreq-1,Nfft).';

U = unpackdof(u,NFreq-1,problem.NInput);
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

States.u     = freq2time(U,NFreq-1,Nfft).';
States.udot  = freq2time(Udot,NFreq-1,Nfft).';
States.uddot = freq2time(Uddot,NFreq-1,Nfft).';

States.w0 = w0;
States.wBase = wBase;

%push through the nl system
o = feval(problem.model,'output',States,hbm,problem).';

%finally convert into a fourier series
O = time2freq(o,NFreq-1,Nfft);

O = packdof(O);