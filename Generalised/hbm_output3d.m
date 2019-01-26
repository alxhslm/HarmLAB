function Op = hbm_output3d(hbm,problem,w0,u,x)
if hbm.options.bUseStandardHBM
    Op = hbm_output(hbm,problem,w0(1),u,x);
    return;
end

NInput = problem.NInput;
NDof   = problem.NDof;

NHarm  = hbm.harm.NHarm;
iSub  = hbm.harm.iSub;
NFreq  = hbm.harm.NFreq;
Nfft   = hbm.harm.Nfft;

iRetain = hbm.harm.iRetain;

%work out the time domain
X = unpackdof(x,NFreq-1,NDof,iRetain);
U = unpackdof(u,NFreq-1,NInput);

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

Op = packdof(O);