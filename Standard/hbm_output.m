function o = hbm_output(hbm,problem,w0,u,x)
NInput = problem.NInput;
NDof = problem.NDof;

NFreq = hbm.harm.NFreq;
Nfft  = hbm.harm.Nfft(1);

iRetain = hbm.harm.iRetain;

if isvector(x) && size(x,1) == hbm.harm.NComp*NDof
    X = unpackdof(x,NFreq-1,NDof,iRetain);
    U = unpackdof(u,NFreq-1,NInput);
else
    X = x;
    U = u;
end

%work out the time domain
States = hbm_states(w0(1),X,U,hbm);

%push through the nl system
o = feval(problem.model,'output',States,hbm,problem).';

%finally convert into a fourier series
O = time2freq(o,NFreq-1,Nfft);

if ndims(x) < 2
    o = packdof(O);
else
    o = O;
end