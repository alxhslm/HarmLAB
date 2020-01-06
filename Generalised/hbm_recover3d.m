function x = hbm_recover3d(hbm,problem,w,u,xnl)
if hbm.options.bUseStandardHBM
    x = hbm_recover(hbm,problem,w,u,xnl);
    return;
end

iRetainNL = hbm.harm.iRetainNL;

r = hbm.harm.rFreqRatio;
w0 = w .* r + hbm.harm.wFreq0;

if ~(isvector(xnl) && size(xnl,1) == hbm.harm.NRetainNL)
    bUnpacked = 1;
    xnl = packdof(xnl,iRetainNL);
    u   = packdof(u);
else
    bUnpacked = 0;
end

if problem.NNL ~= problem.NDof
    f = hbm_nonlinear3d('func',hbm,problem,w0,xnl,u);
    
    A = hbm.lin.Ak + prod(w0)*hbm.lin.Ax;
    B = hbm.lin.Bk + prod(w0)*hbm.lin.Bx;
    for k = 1:2
        A = A + (w0(k)*hbm.lin.Ac{k} + w0(k)^2*hbm.lin.Am{k});
        B = B + (w0(k)*hbm.lin.Bc{k} + w0(k)^2*hbm.lin.Bm{k});
    end
    
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