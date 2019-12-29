function varargout = hbm_balance3d(command,hbm,problem,w,u,x)
if hbm.options.bUseStandardHBM
    [varargout{1:nargout}] = hbm_balance(command,hbm,problem,w,u,x);
    return;
end
NDofTot = hbm.harm.NRetain;

w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

switch command
    case 'func' %F, used by hbm_frf & hbm_bb
        Jx = prod(w0)*hbm.lin.Ax;
        Ju = prod(w0)*hbm.lin.Bx;
        for k = 1:2
            Jx = Jx + (hbm.lin.Ak{k} + w0(k)*hbm.lin.Ac{k} + w0(k)^2*hbm.lin.Am{k});
            Ju = Ju + (hbm.lin.Bk{k} + w0(k)*hbm.lin.Bc{k} + w0(k)^2*hbm.lin.Bm{k});
        end
        cl = Ju*u - Jx*x - hbm.lin.b;
        
        if hbm.bIncludeNL
            cnl = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
        else
            cnl = 0*x;
        end
        c = cl - cnl;
        varargout{1} = c;
    case 'jacob' %dF_dX, used by hbm_frf & hbm_bb
        Jl = -prod(w0)*hbm.lin.Ax;
        for k = 1:2
            Jl = Jl - (hbm.lin.Ak{k} + w0(k)*hbm.lin.Ac{k} + w0(k)^2*hbm.lin.Am{k});
        end
        if hbm.bIncludeNL
            if hbm.options.bAnalyticalDerivs
                [Jx,Jxdot,Jxddot] = hbm_nonlinear3d({'jacobX','jacobXdot','jacobXddot'},hbm,problem,w0,x,u);
                Jnl1 = Jx + w0(1)*Jxdot{1} + w0(2)*Jxdot{2} + w0(1)^2*Jxddot{1} + w0(2)^2*Jxddot{2} + prod(w0)*Jxddot{3};
            else
                c0 = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
                h = 1E-10;
                Jnl2 = zeros(NDofTot);
                x0 = x;
                u0 = u;
                for i = 1:NDofTot
                    x = x0;
                    x(i,:) = x(i,:) + h;
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
        J = Jl - Jnl;
        varargout{1} = J;
    case 'derivW' %dF_dw, used by hbm_frf & hbm_bb
        r = hbm.harm.rFreqRatio;
        Dl = (r(1)*w0(2) + r(2)*w0(1))*(hbm.lin.Bx*u - hbm.lin.Ax*x);
        for k = 1:2
            Dl = Dl + (hbm.lin.Bc{k}*r(k) + 2*r(k)*w0(k)*hbm.lin.Bm{k})*u - (hbm.lin.Ac{k}*r(k) + 2*r(k)*w(k)*hbm.lin.Am{k})*x;
        end
        if hbm.bIncludeNL
            if hbm.dependence.xdot || hbm.dependence.w
                if hbm.options.bAnalyticalDerivs
                    [Jxdot,Jxddot,Judot,Juddot,Dw] = hbm_nonlinear3d({'jacobXdot','jacobXddot','jacobUdot','jacobUddot','derivW'},hbm,problem,w0,x,u);
                    Dnl1 = (r(1)*Jxdot{1} + r(2)*Jxdot{2})*x + (r(1)*Judot{1} + r(2)*Judot{2})*u + r(1)*Dw{1} + r(2)*Dw{2}  + ... 
                            2*w0(1)*Jxddot{1}*x + 2*w0(2)*Jxddot{2}*x + (r(1)*w0(2) + r(2)*w0(1))*Jxddot{3}*x + ...
                            2*w0(1)*Juddot{1}*u + 2*w0(2)*Juddot{2}*u + (r(1)*w0(2) + r(2)*w0(1))*Juddot{3}*u;
                else
                    c0 = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
                    h = 1E-3;
                    for i = 1:2
                        w02 = w0; w02(i) = w02(i) + h;
                        c = hbm_nonlinear3d('func',hbm,problem,w02,x,u);
                        Dw2{i} =  (c-c0)./h;
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
        D = Dl - Dnl;
        varargout{1} = D;
    case 'derivA' %dF_dA, used in hbm_bb
        Dl = 0*x;
        for k = 1:2
            Dl = Dl + (hbm.lin.Bk{k} + w0(k)*hbm.lin.Bc{k} + w0(k)^2*hbm.lin.Bm{k})*u;
        end
        if hbm.bIncludeNL
            [Ju,Judot] = hbm_nonlinear3d({'jacobU','jacobUdot'},hbm,problem,w0,x,u);
            Dnl = (Ju + w0(1)*Judot{1} + w0(2)*Judot{2})*u;
        else
            Dnl = 0*Dl;
        end
        D = Dl - Dnl;
        varargout{1} = D;
    case 'floquet0'
        D0 = -hbm_balance3d('jacob',hbm,problem,w,u,x);
        varargout{1} = D0;
    case 'floquet1' % D2*l^2 + D1*l + D0 = 0, used to solve for floquet multipliers
        D1l = hbm.lin.floquet.D1xdot;
        for k = 1:2
            D1l = D1l + 2*w0(k)*hbm.lin.floquet.D1xddot{k};
        end
        if hbm.bIncludeNL
            [D1xdot,D1xddot] = hbm_nonlinear3d({'floquet1xdot','floquet1xddot'},hbm,problem,w0,x,u);
            D1nl = D1xdot;
            for k = 1:2
                D1nl = D1nl + 2*w0(k)*D1xddot{k};
            end
        else
            D1nl = 0*D1l;
        end
        varargout{1} = D1l - D1nl;
    case 'floquet2'
        D2l = hbm.lin.floquet.D2;
        if hbm.bIncludeNL
            D2nl = hbm_nonlinear3d({'floquet2'},hbm,problem,w0,x,u);
        else
            D2nl = 0*D2l;
        end
        varargout{1} = D2l - D2nl;
        %% All cases below are currently unused - they are required for computing the resonance
    case 'jacobW' %d2F_dwdX, to be used in hbm_bb
        Dl = - (hbm.lin.Ac + 2*w0*hbm.lin.Am);
        if hbm.bIncludeNL
            Jxdot = hbm_nonlinear3d('jacobXdot',hbm,problem,w0,x,u);
            Dnl = Jxdot;
        else
            Dnl = 0*Dl;
        end
        D = Dnl + Dl;
        varargout{1} = D;
    case 'derivW2' %d2F_dw2, to be used in hbm_bb
        Dl = 2*hbm.lin.Bm*u - 2*hbm.lin.Am*x;
        Dnl = 0*Dl;
        D = Dnl + Dl;
        varargout{1} = D;
    case 'derivWA' %d2F_dwdA, to be used in hbm_bb
        Dl = (hbm.lin.Bc + 2*w0*hbm.lin.Bm)*u;
        if hbm.bIncludeNL
            Judot = hbm_nonlinear3d({'jacobUdot'},hbm,problem,w0,x,u);
            Dnl = Judot*u;
        else
            Dnl = 0*Dl;
        end
        D = Dnl + Dl;
        varargout{1} = D;
end