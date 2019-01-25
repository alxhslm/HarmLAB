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
rBase  = hbm.harm.rFreqBase;

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;

%unpack the inputs
w = kHarm*(w0.*rBase)';
t1 = (0:Nfft(1)-1)/Nfft(1)*2*pi/(rBase(1)*w0(1));
t2 = (0:Nfft(2)-1)/Nfft(2)*2*pi/(rBase(2)*w0(2));
[t1,t2] = ndgrid(t1,t2);
t = [t1(:) t2(:)];

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
        x     = freq2time3d(X,NHarm,hbm.harm.iSub,Nfft);
        xdot  = freq2time3d(Xdot,NHarm,hbm.harm.iSub,Nfft);
        xddot = freq2time3d(Xddot,NHarm,hbm.harm.iSub,Nfft);
        
        %create the vector of inputs
        u     = freq2time3d(U,NHarm,hbm.harm.iSub,Nfft);
        udot  = freq2time3d(Udot,NHarm,hbm.harm.iSub,Nfft);
        uddot = freq2time3d(Uddot,NHarm,hbm.harm.iSub,Nfft);

    case 'mat'       
        %create the time series from the fourier series
        x     = real(hbm.nonlin.IFFT*X);
        xdot  = real(hbm.nonlin.IFFT*Xdot);
        xddot = real(hbm.nonlin.IFFT*Xddot);
        
        %create the vector of inputs
        u     = real(hbm.nonlin.IFFT*U);
        udot  = real(hbm.nonlin.IFFT*Udot);
        uddot = real(hbm.nonlin.IFFT*Uddot);
end

%push through the nl system
o = feval(problem.model,'output',t', x.',xdot.',xddot.',u.',udot.',uddot.',hbm,problem,w0).';

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        %put back into hypertime
        O = time2freq3d(o,NHarm,hbm.harm.iSub,Nfft);
    case 'mat'
        %finally convert into a fourier series
        O = hbm.nonlin.FFT*o;
end

Op = packdof(O);