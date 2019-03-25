function States = hbm_states3d(w0,X,U,hbm) 
iSub  = hbm.harm.iSub;
NHarm  = hbm.harm.NHarm;
Nfft  = hbm.harm.Nfft;
kHarm = hbm.harm.kHarm;
wBase = hbm.harm.rFreqBase.*w0;

%unpack the inputs
w = kHarm*wBase';
t1 = (0:Nfft(1)-1)/Nfft(1)*2*pi/wBase(1);
t2 = (0:Nfft(2)-1)/Nfft(2)*2*pi/wBase(2);
[t1,t2] = ndgrid(t1,t2);
States.t = [t1(:) t2(:)].';

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
switch hbm.options.aft_method
    case 'fft'
        %create the time series from the fourier series
        States.x     = freq2time3d(X,NHarm,iSub,Nfft).';
        States.xdot  = freq2time3d(Xdot,NHarm,iSub,Nfft).';
        States.xddot = freq2time3d(Xddot,NHarm,iSub,Nfft).';
        
        %create the vector of inputs
        States.u     = freq2time3d(U,NHarm,iSub,Nfft).';
        States.udot  = freq2time3d(Udot,NHarm,iSub,Nfft).';
        States.uddot = freq2time3d(Uddot,NHarm,iSub,Nfft).';
        
    case 'mat'
        %create the time series from the fourier series
        States.x     = real(hbm.nonlin.IFFT*X).';
        States.xdot  = real(hbm.nonlin.IFFT*Xdot).';
        States.xddot = real(hbm.nonlin.IFFT*Xddot).';
        
        %create the vector of inputs
        States.u     = real(hbm.nonlin.IFFT*U).';
        States.udot  = real(hbm.nonlin.IFFT*Udot).';
        States.uddot = real(hbm.nonlin.IFFT*Uddot).';
end