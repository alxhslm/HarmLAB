function x = hbm_recover(hbm,problem,w,u,xnl)
NNL = problem.NNL;
iRetainNL = hbm.harm.iRetainNL;

ii = find(hbm.harm.NHarm ~= 0);
r = hbm.harm.rFreqRatio(ii);
w0 = w * r + hbm.harm.wFreq0(ii);

if ~(isvector(xnl) && size(xnl,1) == hbm.harm.NComp*NNL)
    bUnpacked = 1;
    xnl = packdof(xnl,iRetainNL);
    u   = packdof(u);
else
    bUnpacked = 0;
end

if problem.NNL ~= problem.NDof
    f = hbm_nonlinear('func',hbm,problem,w0,xnl,u);
    
    A = (hbm.lin.Ak{ii} + w0*hbm.lin.Ac{ii} + w0^2*hbm.lin.Am{ii});
    B = (hbm.lin.Bk{ii} + w0*hbm.lin.Bc{ii} + w0^2*hbm.lin.Bm{ii});
    
    
    iLin = hbm.harm.iLin;
    iNL = hbm.harm.iNL;

    App = A(iLin,iLin);
    Apq = A(iLin,iNL);
    
    x = zeros(hbm.harm.NRetain,1);
    x(iLin) = App\(B(iLin,:)*u - hbm.lin.b(iLin) - f(iLin) - Apq*xnl);
    x(iNL) = xnl;
else
    x = xnl;
end

if bUnpacked
    x = unpackdof(x,hbm.harm.NHarm,problem.NDof,hbm.harm.iRetain);
end