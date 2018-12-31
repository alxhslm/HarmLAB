function varargout = hbm_objective(part,hbm,problem,w0,x,u)

NDof = problem.NDof;
if problem.res.iHarm > 1
    iRe = problem.res.iDof + 2*(problem.res.iHarm - 2)*NDof + NDof;
    iIm = problem.res.iDof + 2*(problem.res.iHarm - 2)*NDof + 2*NDof;
else
    iRe = problem.res.iDof;
    iIm = [];
end

W = w0(1)/hbm.harm.rFreqRatio(1);
r = hbm.harm.kHarm*(hbm.harm.rFreqBase.*hbm.harm.rFreqRatio)';
D = [0 -1;
    1  0];
D = blkdiag(0,kron(diag(r(2:end)),D));

%% output
Wx = kron(D,eye(problem.NDof));
switch problem.res.output
    case 'none'
        Fnl = 0*x;
    case 'fnl'
        Fnl = hbm_nonlinear3d({'func'},hbm,problem,w0,x,u);
    case 'x'
        Fnl = x;
    case 'xdot'
        Fnl = W*Wx*x;
    case 'xddot'
        Fnl = W^2*Wx*Wx*x;
end

if ~isempty(iIm)
    Fb = sqrt(Fnl(iRe)^2 + Fnl(iIm).^2);
else
    Fb = Fnl(iRe);
end

%% input
Wu = kron(D,eye(problem.NInput));
switch problem.res.input
    case 'unity'
        Fex = 0*u+1;
    case 'fe'
        Ju = prod(w0)*hbm.lin.Bx;
        for k = 1:2
            Ju = Ju + (hbm.lin.Bk{k} + w0(k)*hbm.lin.Bc{k} + w0(k)^2*hbm.lin.Bm{k});
        end
        Fex = Ju*u -  hbm.lin.b;
    case 'u'
        Fex = u;
    case 'xdot'
        Fex = W * Wu*u;
    case 'xddot'
        Fex = W^2 * Wu*Wu*u;
end

if ~isempty(iIm)
    Fe = sqrt(Fex(iRe)^2 + Fex(iIm).^2);
else
    Fe = Fex(iRe);
end

if ~iscell(part)
    part = {part};
end
varargout = cell(1,length(part));

for i = 1:length(part)
    switch part{i}
        case 'func'
            H = Fb./Fe;
            varargout{i} = H;
        case 'jacobX'
            
            %nl
            switch problem.res.output
                case 'none'
                    Jx = 0*Wx;
                case 'fnl'
                    Jx = hbm_nonlinear3d({'jacobX'},hbm,problem,w0,x,u);
                case 'x'
                    Jx = eye(problem.NDof*hbm.harm.NComp);
                case 'xdot'
                    Jx = Wx;
                case 'xddot'
                    Jx = Wx*Wx;
            end
            
            if ~isempty(iIm)
                dFbdx = (Fnl(iRe)*Jx(iRe,:) + Fnl(iIm)*Jx(iIm,:))./(Fb + eps);
            else
                dFbdx = (Fnl(iRe)*Jx(iRe,:))./(Fb + eps);
            end
            
            %excitation
            dFedx = 0;
            
            %put it all together
            dHdx = (dFbdx.*Fe - Fb.*dFedx)./Fe.^2;
            
            varargout{i} = dHdx;
        case 'derivW'
            r = hbm.harm.rFreqRatio;
            
            %nl
            switch problem.res.output
                case 'none'
                    Dw_nl = 0*x;
                case 'fnl'
                    Dw = hbm_nonlinear3d({'derivW'},hbm,problem,w0,x,u);
                    if ~iscell(Dw)
                        Dw = {Dw,0*Dw};
                    end
                    Dw_nl = 0*Dw{1};
                    for k = 1:2
                        Dw_nl = Dw_nl + r(k)*Dw{k};
                    end
                case 'x'
                    Dw_nl = zeros(problem.NDof*hbm.harm.NComp,1);
                case 'xdot'
                    Dw_nl = (Wx*x);
                case 'xddot'
                    Dw_nl = 2*W*(Wx*Wx*x);
            end
            if ~isempty(iIm)
                dFbdw = (Fnl(iRe)*Dw_nl(iRe) + Fnl(iIm)*Dw_nl(iIm))./(Fb + eps);
            else
                dFbdw = (Fnl(iRe)*Dw_nl(iRe))./(Fb + eps);
            end
            
            %excitation
            switch problem.res.input
                case 'unity'
                    Dw_u = 0*u;
                case 'fe'
                    Dw_u = 2*prod(r)*W*hbm.lin.Bx*u;
                    for k = 1:2
                        Dw_u = Dw_u + (hbm.lin.Bc{k}*r(k) + 2*r(k)^2*W*hbm.lin.Bm{k})*u;
                    end
                case 'u'
                    Dw_u = zeros(problem.NInput*hbm.harm.NComp,1);
                case 'udot'
                    Dw_u = (Wu*u);
                case 'uddot'
                    Dw_u = 2*W*(Wu*Wu*u);
            end
            if ~isempty(iIm)
                dFedw = (Fex(iRe)*Dw_u(iRe) + Fex(iIm)*Dw_u(iIm))./Fe;
            else
                dFedw = (Fex(iRe)*Dw_u(iRe))./Fe;
            end
            
            dHdw = (dFbdw.*Fe - Fb.*dFedw)./Fe.^2;
            
            varargout{i} = dHdw;
    end
end