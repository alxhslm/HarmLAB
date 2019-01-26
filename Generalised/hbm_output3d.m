function Op = hbm_output3d(hbm,problem,w0,u,x)
if hbm.options.bUseStandardHBM
    Op = hbm_output(hbm,problem,w0(1),u,x);
    return;
end

NInput = problem.NInput;
NDof   = problem.NDof;

NHarm  = hbm.harm.NHarm;
NComp  = hbm.harm.NComp;
NFreq  = hbm.harm.NFreq;
Nfft   = hbm.harm.Nfft;
kHarm  = hbm.harm.kHarm;
wBase  = hbm.harm.rFreqBase.*w0;

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;

%unpack the inputs
w = kHarm*wBase';
t1 = (0:Nfft(1)-1)/Nfft(1)*2*pi/wBase(1);
t2 = (0:Nfft(2)-1)/Nfft(2)*2*pi/wBase(2);
[t1,t2] = ndgrid(t1,t2);
States.t = [t1(:) t2(:)].';

%now compute the fourier coefficients for xdot, xddot
X = unpackdof(x,NFreq-1,NDof,iRetain);
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

U = unpackdof(u,NFreq-1,problem.NInput);
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

switch hbm.options.aft_method
    case 'fft'        
        %create the time series from the fourier series
        States.x     = freq2time3d(X,NHarm,hbm.harm.iSub,Nfft)';
        States.xdot  = freq2time3d(Xdot,NHarm,hbm.harm.iSub,Nfft).';
        States.xddot = freq2time3d(Xddot,NHarm,hbm.harm.iSub,Nfft)';
        
        %create the vector of inputs
        States.u     = freq2time3d(U,NHarm,hbm.harm.iSub,Nfft)';
        States.udot  = freq2time3d(Udot,NHarm,hbm.harm.iSub,Nfft)';
        States.uddot = freq2time3d(Uddot,NHarm,hbm.harm.iSub,Nfft)';

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

States.w0 = w0;
States.wBase = wBase;

%push through the nl system
o = feval(problem.model,'output',States,hbm,problem);

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        %put back into hypertime
        O = time2freq3d(o.',NHarm,hbm.harm.iSub,Nfft);
    case 'mat'
        %finally convert into a fourier series
        O = hbm.nonlin.FFT*o.';
end

Op = packdof(O);