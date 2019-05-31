function x = hbm_recover(hbm,problem,w0,u,xnl)
NNL = problem.NNL;
iRetainNL = hbm.harm.iRetainNL;

if ~(isvector(xnl) && size(xnl,1) == hbm.harm.NComp*NNL)
    bUnpacked = 1;
    xnl = packdof(xnl,iRetainNL);
    u   = packdof(u);
else
    bUnpacked = 0;
end

if problem.NNL ~= problem.NDof
    f = hbm_nonlinear('func',hbm,problem,w0,xnl,u);
    
    A = (hbm.lin.Ak{1} + w0*hbm.lin.Ac{1} + w0^2*hbm.lin.Am{1});
    B  = (hbm.lin.Bk{1} + w0*hbm.lin.Bc{1} + w0^2*hbm.lin.Bm{1});
    
    
    iLin = hbm.harm.iLin;
    iNL = hbm.harm.iNL;
    
    App = A(iLin,iLin);
    Apq = A(iLin,iNL);
    
    x = zeros(hbm.harm.NRetain,1);
    x(iLin) = App\(B(iLin,:)*u - f(iLin) - Apq*xnl);
    x(iNL) = xnl;
end

if bUnpacked
    x = unpackdof(x,hbm.harm.NHarm,problem.NDof,hbm.harm.iRetain);
end