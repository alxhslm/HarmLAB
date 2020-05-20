function o = hbm_output(hbm,problem,w,u,x)
NInput = problem.NInput;
NDof = problem.NDof;
NHarm  = hbm.harm.NHarm;

ii = find(NHarm ~= 0);
r = hbm.harm.rFreqRatio(ii);
w0 = w * r + hbm.harm.wFreq0(ii);

if isvector(x) && size(x,1) == hbm.harm.NComp*problem.NDof
    X = unpackdof(x,NFreq-1,NDof);
    U = unpackdof(u,NFreq-1,NInput);
else
    X = x;
    U = u;
end

%work out the time domain
States = hbm_states(w0,X,U,hbm);

%push through the nl system
o = feval(problem.model,'output',States,hbm,problem);

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        %put back into hypertime
        O = time2freq(o.',NHarm,NFreq-1,Nfft);
    case 'mat'
        %finally convert into a fourier series
        O = hbm.nonlin.FFT*o.';
end

if ndims(x) < 2
    o = packdof(O);
else
    o = O;
end