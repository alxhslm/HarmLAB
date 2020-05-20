function States = hbm_states(w0,X,U,hbm)
ii = find(hbm.harm.NHarm ~= 0);

NFreq = hbm.harm.NFreq;
Nfft  = hbm.harm.Nfft(ii);
kHarm = hbm.harm.kHarm(:,ii);
wBase = hbm.harm.rFreqBase(ii)*w0;

%unpack the inputs
w = kHarm*wBase;
States.t = (0:Nfft-1)/Nfft*2*pi/wBase;

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%precompute the external inputs
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

States.w0 = w0;
States.wBase = wBase;

%create the time series from the fourier series
States.x     = real(hbm.nonlin.IFFT*X).';
States.xdot  = real(hbm.nonlin.IFFT*Xdot).';
States.xddot = real(hbm.nonlin.IFFT*Xddot).';

%create the vector of inputs
States.u     = real(hbm.nonlin.IFFT*U).';
States.udot  = real(hbm.nonlin.IFFT*Udot).';
States.uddot = real(hbm.nonlin.IFFT*Uddot).';
