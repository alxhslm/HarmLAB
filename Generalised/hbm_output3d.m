function o = hbm_output3d(hbm,problem,w,u,x)
if hbm.options.bUseStandardHBM
    o = hbm_output(hbm,problem,w,u,x);
    return;
end

NInput = problem.NInput;
NDof   = problem.NDof;

r = hbm.harm.rFreqRatio;
w0 = w .* r + hbm.harm.wFreq0;

if size(x,1) == hbm.harm.NFreq && size(x,2) == problem.NDof
    X = x;
    U = u;
elseif  isvector(x) && size(x,1) == hbm.harm.NComp*problem.NDof
    X = unpackdof(x,NFreq-1,NDof);
    U = unpackdof(u,NFreq-1,NInput);
end

%work out the time domain
States = hbm_states3d(w0,X,U,hbm);

%push through the nl system
o = feval(problem.model,'output',States,hbm,problem);

%finally convert into a fourier series
O = hbm.nonlin.FFT*o.';

if size(x,1) == hbm.harm.NFreq && size(x,2) == problem.NDof
    o = O;
else
    o = packdof(O);    
end