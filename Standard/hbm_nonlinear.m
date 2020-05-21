function varargout = hbm_nonlinear(command,hbm,problem,w0,Xp,Up)
NInput = problem.NInput;
NDof = problem.NDof;
NNL  = problem.NNL;

ii = find(hbm.harm.NHarm ~= 0);

NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
NDofTot = NComp*NDof;
NNLTot  = NComp*NNL;
NInputTot = NComp*NInput;

%work out the time domain
X = zeros(NFreq,NDof);
X(:,problem.iNL) = unpackdof(Xp,NFreq-1,NNL);
U = unpackdof(Up,NFreq-1,NInput);

%get the time series
States = hbm_states(w0,X,U,hbm);

%push through the nl system
States.f = feval(problem.model,'nl' ,States,hbm,problem);

%finally convert into the time domain
F = hbm.nonlin.FFT*States.f.';
Fp = packdof(F);

if ~iscell(command)
    command = {command};
end

varargout = cell(1,length(command));

ijacobx   = hbm.nonlin.hbm.ijacobx;
ijacobxnl = hbm.nonlin.hbm.ijacobxnl;
ijacobu   = hbm.nonlin.hbm.ijacobu;

for o = 1:length(command)
    switch command{o}
        case 'func'
            %all done
            varargout{o} = Fp;
        case 'derivW' %df_dw
            if ~hbm.dependence.w
                Dw = zeros(NDofTot,1);
            else
                dfhbm_dw = hbm_derivatives('nl' ,'w',States,hbm,problem);
                if isfield(problem,'derivW')
                    Dw = feval(problem.derivW,dfhbm_dw,hbm,problem);
                else
                    Dw = packdof(hbm.nonlin.FFT*dfhbm_dw{ii});
                end
            end
            varargout{o} = Dw;
        case 'jacobX'  %df_dX = Jx*X
            if ~hbm.dependence.x
                Jx = zeros(NDofTot,NNLTot);
            else
                States.df_dx = hbm_derivatives('nl','x',States,hbm,problem);
                
                if isfield(problem,'jacobX')
                    Jx = feval(problem.jacobX,States,hbm,problem);
                else
                    Jx = sum(hbm.nonlin.hbm.Jx.*States.df_dx(ijacobx,ijacobxnl,:),3);
                end
            end
            varargout{o} = Jx;
        case 'jacobU' %df_dU = Ju*U
            if ~hbm.dependence.u
                Ju = zeros(NDofTot,NInputTot);
            else
                States.df_du = hbm_derivatives('nl' ,'u',States,hbm,problem);
                
                if isfield(problem,'jacobU')
                    Ju = feval(problem.jacobU,States,problem);
                else
                    Ju = sum(hbm.nonlin.hbm.Ju.*States.df_du(ijacobx,ijacobu,:),3);
                end
            end
            varargout{o} = Ju;
        case 'jacobXdot' %df_dX = w*Jxdot*X
             if ~hbm.dependence.xdot
                 Jxdot = zeros(DofTot,NNLTot);
             else
                 States.df_dxdot = hbm_derivatives('nl' ,'xdot',States,hbm,problem);
                 
                 if isfield(problem,'jacobXdot')
                     Jxdot = feval(problem.jacobXdot,States,hbm,problem);
                 else
                     Jxdot = sum(hbm.nonlin.hbm.Jxdot{ii}.*States.df_dxdot(ijacobx,ijacobxnl,:),3);
                 end
             end
            varargout{o} = Jxdot;
        case 'jacobUdot' %df_dU = w*Judot*U
            if ~hbm.dependence.udot
                Judot = zeros(NDofTot,NInputTot);
            else
                States.df_dudot = hbm_derivatives('nl' ,'udot',States,hbm,problem);
                
                if isfield(problem,'jacobUdot')
                    Judot = feval(problem.jacobUdot,States,hbm,problem);
                else
                    Judot = sum(hbm.nonlin.hbm.Judot{ii}.*States.df_dudot(ijacobx,ijacobu,:),3);
                end
            end
            varargout{o} = Judot;
         case 'jacobXddot' %df_dX = w^2*Jxddot*X
            if ~hbm.dependence.xddot
                Jxddot = zeros(NDofTot,NNLTot);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'jacobXddot')
                    Jxddot = feval(problem.jacobXddot,States,hbm,problem);
                else
                    Jxddot = sum(hbm.nonlin.hbm.Jxddot{ii}.*States.df_dxddot(ijacobx,ijacobxnl,:),3);
                end
            end
            varargout{o} = Jxddot;
        case 'jacobUddot' %df_dU = w^2*Juddot*U
            if ~hbm.dependence.uddot
                Juddot = zeros(NDofTot,NInputTot);
            else
                States.df_duddot = hbm_derivatives('nl' ,'uddot',States,hbm,problem);
                
                if isfield(problem,'jacobUddot')
                    Juddot = feval(problem.jacobUddot,States,hbm,problem);
                else
                    Juddot = sum(hbm.nonlin.hbm.Juddot{ii}.*States.df_duddot(ijacobx,ijacobu,:),3);
                end
            end
            varargout{o} = Juddot;
       case 'floquet1xdot'
           if ~hbm.dependence.xdot
                D1 = zeros(NDofTot,NNLTot);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);
                
                if isfield(problem,'floquet1xdot')
                    D1 = feval(problem.floquet1xdot,States,hbm,problem);
                else
                    D1 = sum(hbm.nonlin.hbm.Jx.*States.df_dxdot(ijacobx,ijacobxnl,:),3);
                end
            end
            varargout{o} = D1;
        case 'floquet1xddot'
           if ~hbm.dependence.xddot
                D1dd = zeros(NDofTot,NNLTot);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'floquet1xddot')
                    D1dd = feval(problem.floquet1xddot,States,hbm,problem);
                else
                    D1dd = sum(hbm.nonlin.hbm.Jxdot{ii}.*States.df_dxddot(ijacobx,ijacobxnl,:),3);
                end
            end
            varargout{o} = D1dd; 
       case 'floquet2'
           if ~hbm.dependence.xddot
                D2 = zeros(NDofTot,NNLTot);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'floquet2')
                    D2 = feval(problem.floquet2,States,hbm,problem);
                else
                    D2 = sum(hbm.nonlin.hbm.Jx.*States.df_dxddot(ijacobx,ijacobxnl,:),3);
                end
            end
            varargout{o} = D2;
    end
end