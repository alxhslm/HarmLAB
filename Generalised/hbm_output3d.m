function o = hbm_output3d(hbm,problem,w,u,x)
if hbm.options.bUseStandardHBM
    o = hbm_output(hbm,problem,w,u,x);
    return;
end

NInput = problem.NInput;
NDof   = problem.NDof;

NHarm  = hbm.harm.NHarm;
iSub  = hbm.harm.iSub;
NFreq  = hbm.harm.NFreq;
Nfft   = hbm.harm.Nfft;

iRetain = hbm.harm.iRetain;

r = hbm.harm.rFreqRatio;
w0 = w .* r + hbm.harm.wFreq0;

if isvector(x)
    X = unpackdof(x,NFreq-1,NDof,iRetain);
    U = unpackdof(u,NFreq-1,NInput);
else
    X = x;
    U = u;
end

%work out the time domain
States = hbm_states3d(w0,X,U,hbm);

%push through the nl system
o = feval(problem.model,'output',States,hbm,problem);

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        %put back into hypertime
        O = time2freq3d(o.',NHarm,iSub,Nfft);
    case 'mat'
        %finally convert into a fourier series
        O = hbm.nonlin.FFT*o.';
end

if isvector(x)
    o = packdof(O);
else
    o = O;
end