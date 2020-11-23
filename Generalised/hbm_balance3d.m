function varargout = hbm_balance3d(command,hbm,problem,w,u,x)
if hbm.options.bUseStandardHBM
    [varargout{1:nargout}] = hbm_balance(command,hbm,problem,w,u,x);
    return;
end
NDofTot = hbm.harm.NComp*problem.NDof;

r = hbm.harm.rFreqRatio;
w0 = w .* r + hbm.harm.wFreq0;

A = hbm.lin.Ak + prod(w0)*hbm.lin.Ax;
B = hbm.lin.Bk + prod(w0)*hbm.lin.Bx;
dAdw = (r(1)*w0(2) + r(2)*w0(1))*hbm.lin.Ax;
dBdw = (r(1)*w0(2) + r(2)*w0(1))*hbm.lin.Bx;
for k = 1:2
    A = A + (w0(k)*hbm.lin.Ac{k} + w0(k)^2*hbm.lin.Am{k});
    B = B + (w0(k)*hbm.lin.Bc{k} + w0(k)^2*hbm.lin.Bm{k});

    dAdw  = dAdw  + r(k)*(hbm.lin.Ac{k} + 2*w0(k)*hbm.lin.Am{k});
    dBdw  = dBdw  + r(k)*(hbm.lin.Bc{k} + 2*w0(k)*hbm.lin.Bm{k});
end

switch command
    case 'func' %F, used by hbm_frf & hbm_bb
        cl = B*u - hbm.lin.b - A*x;
        if hbm.bIncludeNL
            cnl = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
        else
            cnl = 0*cl;
        end
        c = (cl - cnl);
        varargout{1} = c;
    case 'jacob' %dF_dX, used by hbm_frf & hbm_bb
        Jl = -A;
        if hbm.bIncludeNL
            if hbm.options.bAnalyticalDerivs
                [Jx,Jxdot,Jxddot] = hbm_nonlinear3d({'jacobX','jacobXdot','jacobXddot'},hbm,problem,w0,x,u);
                Jnl1 = Jx + w0(1)*Jxdot{1} + w0(2)*Jxdot{2} + w0(1)^2*Jxddot{1} + w0(2)^2*Jxddot{2} + prod(w0)*Jxddot{3};
            else
                c0 = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
                h = 1E-12;
                Jnl2 = zeros(NDofTot,NDofTot);
                x0 = x;
                u0 = u;
                for i = 1:NDofTot
                    x = x0;
                    x(i) = x(i) + h;
                    c = hbm_nonlinear3d('func',hbm,problem,w0,x,u0);
                    Jnl2(:,i) = (c-c0)./h;
                end
            end
            if hbm.options.bAnalyticalDerivs
                Jnl = Jnl1;
            else
                Jnl = Jnl2;
            end
        else
            Jnl = 0*Jl;
        end
        J = (Jl - Jnl);
        varargout{1} = J;
    case 'derivW' %dF_dw, used by hbm_frf & hbm_bb
        Dl = dBdw*u - dAdw*x;
        cl = B*u - hbm.lin.b - A*x;
        if hbm.bIncludeNL
            if hbm.dependence.xdot || hbm.dependence.w
                if hbm.options.bAnalyticalDerivs
                    [Jxdot,Jxddot,Judot,Juddot,Dw] = hbm_nonlinear3d({'jacobXdot','jacobXddot','jacobUdot','jacobUddot','derivW'},hbm,problem,w0,x,u);
                    Dxdot  = (r(1)*Jxdot{1} + r(2)*Jxdot{2})*x;
                    Dxddot = (2*r(1)*w0(1)*Jxddot{1} + 2*r(2)*w0(2)*Jxddot{2} + (r(1)*w0(2) + r(2)*w0(1))*Jxddot{3})*x;
                    Dudot  = (r(1)*Judot{1} + r(2)*Judot{2})*u;
                    Duddot = (2*r(1)*w0(1)*Juddot{1} + 2*r(2)*w0(2)*Juddot{2} + (r(1)*w0(2) + r(2)*w0(1))*Juddot{3})*u;
                    Dw = r(1)*Dw{1} + r(2)*Dw{2};
                    Dnl1 = Dxdot + Dxddot + Dudot + Duddot + Dw;
                else
                    cnl = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
                    h = 1E-10;
                    for i = 1:2
                        w02 = w0; w02(i) = w02(i) + h;
                        c = hbm_nonlinear3d('func',hbm,problem,w02,x,u);
                        Dw2{i} =  (c-cnl)./h;
                    end
                    Dnl2 = Dw2{1}*r(1) + Dw2{2}*r(2);
                end
                if hbm.options.bAnalyticalDerivs
                    Dnl = Dnl1;
                else
                    Dnl = Dnl2;
                end
            else
                Dnl = 0*Dl;
            end
        else
            Dnl = 0*Dl;
        end
        D = (Dl - Dnl);
        varargout{1} = D;
    case 'derivA' %dF_dA, used in hbm_bb
        Dl = B*u;
        if hbm.bIncludeNL
            [Ju,Judot,Juddot] = hbm_nonlinear3d({'jacobU','jacobUdot','jacobUddot'},hbm,problem,w0,x,u);
            Dnl = (Ju + w0(1)*Judot{1} + w0(2)*Judot{2} +w0(1)^2*Juddot{1} + w0(2)^2*Juddot{2} + prod(w0)*Juddot{3})*u;
        else
            Dnl = 0*Dl;
        end
        D = (Dl - Dnl);
        varargout{1} = D;
    case 'floquet0'
        D0 = -hbm_balance3d('jacob',hbm,problem,w,u,x);
        varargout{1} = D0;
    case 'floquet1'
        D1l = hbm.lin.floquet.D1xdot + hbm.lin.floquet.D1Gxdot*w0(1) + 2*(w0(1)*hbm.lin.floquet.D1xddot{1} + w0(2)*hbm.lin.floquet.D1xddot{2});
        if hbm.bIncludeNL
            [D1xdot,D1xddot] = hbm_nonlinear3d({'floquet1xdot','floquet1xddot'},hbm,problem,w0,x,u);
            D1nl = D1xdot + 2*w0(1)*D1xddot{1} + 2*w0(2)*D1xddot{2};
        else
            D1nl = 0*D1l;
        end
        varargout{1} = (D1l - D1nl);
    case 'floquet2'
        D2l = hbm.lin.floquet.D2;
        if hbm.bIncludeNL
            D2nl = hbm_nonlinear3d({'floquet2'},hbm,problem,w0,x,u);
        else
            D2nl = 0*D2l;
        end
        varargout{1} = (D2l - D2nl);
end