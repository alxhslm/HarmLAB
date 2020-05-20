function varargout = hbm_nonlinear3d(command,hbm,problem,w0,Xp,Up)
if hbm.options.bUseStandardHBM
    [varargout{1:nargout}] = hbm_nonlinear(command,hbm,problem,w0(1),Xp,Up);
    return;
end

NInput = problem.NInput;
NNL = problem.NNL;
NDof = problem.NDof;

NComp  = hbm.harm.NComp;
NFreq  = hbm.harm.NFreq;
Nfft   = hbm.harm.Nfft;
kHarm  = hbm.harm.kHarm;

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;
iRetainNL = hbm.harm.iRetainNL;
NRetainNL = hbm.harm.NRetainNL;
NDofTot = NComp*NDof;
NNLTot  = NComp*NNL;
NInputTot = NComp*NInput;
iNL = problem.iNL;

%work out the time domain
X = zeros(NFreq,NDof);
X(:,iNL) = unpackdof(Xp,NFreq-1,NNL,iRetainNL);

U = unpackdof(Up,NFreq-1,NInput);

%get the time series
States = hbm_states3d(w0,X,U,hbm);

%push through the nl system
States.f = feval(problem.model,'nl',States,hbm,problem);

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        %put back into hypertime
        F = time2freq3d(States.f.',NHarm,hbm.harm.iSub,Nfft);
    case 'mat'
        %finally convert into a fourier series
        F = hbm.nonlin.FFT*States.f.';
end

Fp = packdof(F,iRetain);

if ~iscell(command)
    command = {command};
end

varargout = cell(1,length(command));

ijacobx = hbm.nonlin.hbm.ijacobx;
ijacobxnl = hbm.nonlin.hbm.ijacobxnl;
ijacobu = hbm.nonlin.hbm.ijacobu;

for o = 1:length(command)
    switch command{o}
        case 'func'
            %all done
            varargout{o} = Fp;
        case 'derivW' %df_dw
            if ~hbm.dependence.w
                Dw = repmat({zeros(NRetain,1)},1,2);
            else
                dfhbm_dw = hbm_derivatives('nl','w',States,hbm,problem);
                if isfield(problem,'derivW')
                    Dw = feval(problem.derivW,dfhbm_dw,hbm,problem);
                else
                    for n = 1:2
                        Dw{n} = packdof(hbm.nonlin.FFT*dfhbm_dw{n},iRetain);
                    end
                end
            end
            varargout{o} = Dw;
        case 'jacobX'  %df_dX = Jx*X
            if ~hbm.dependence.x
                Jx = zeros(NRetain,NRetainNL);
            else
                States.df_dx = hbm_derivatives('nl','x',States,hbm,problem);
                
                if isfield(problem,'jacobX')
                    Jx = feval(problem.jacobX,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jx = sum(hbm.nonlin.hbm.Jx.*States.df_dx(ijacobx,ijacobxnl,:),3);
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Jx = zeros(NDofTot,NNLTot);
                            Jx(1:NDof,1:NNL) = mean(States.df_dx(:,iNL,:),3);
                            for l = 1:(NFreq-1)
                                [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                Jx(1:NDof,(2*l-1)*NNL+(1:NNL)) =  sum( (States.df_dx(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                                Jx(1:NDof,(2*l-0)*NNL+(1:NNL)) =  sum(-(States.df_dx(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                Jx((2*k-1)*NDof+(1:NDof),1:NNL) =  2*sum(States.df_dx(:,iNL,:).*ck/(Nfft(1)*Nfft(2)),3);
                                Jx((2*k-0)*NDof+(1:NDof),1:NNL) = -2*sum(States.df_dx(:,iNL,:).*sk/(Nfft(1)*Nfft(2)),3);
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) =  2*sum( (States.df_dx(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) =  2*sum(-(States.df_dx(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) = -2*sum( (States.df_dx(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) = -2*sum(-(States.df_dx(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                            Jx = Jx(iRetain,iRetainNL);
                    end
                end
            end
            varargout{o} = Jx;
        case 'jacobU' %df_dU = Ju*U
            if ~hbm.dependence.u
                Ju = zeros(NRetain,NInputTot);
            else
                States.df_du = hbm_derivatives('nl','u',States,hbm,problem);
                
                if isfield(problem,'jacobU')
                    Ju = feval(problem.jacobU,States.df_du,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Ju = sum(hbm.nonlin.hbm.Ju.*States.df_du(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Ju = zeros(NDofTot,NInputTot);
                            Ju(1:NDof,1:NInput) = mean(States.df_du,3);
                            for l = 1:(NFreq-1)
                                [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                Ju(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum( (States.df_du.*cl)/(Nfft(1)*Nfft(2)),3);
                                Ju(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum(-(States.df_du.*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                Ju((2*k-1)*NDof+(1:NDof),1:NInput) =  2*sum(States.df_du.*ck/(Nfft(1)*Nfft(2)),3);
                                Ju((2*k-0)*NDof+(1:NDof),1:NInput) = -2*sum(States.df_du.*sk/(Nfft(1)*Nfft(2)),3);
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Ju((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum( (States.df_du.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Ju((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum(-(States.df_du.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Ju((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum( (States.df_du.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    Ju((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum(-(States.df_du.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                            Ju = Ju(iRetain,:);
                    end
                end
            end
            varargout{o} = Ju;
        case 'jacobXdot' %df_dX = w*Jxdot*X
            if ~hbm.dependence.xdot
                Jxdot = repmat({zeros(NRetain,NRetainNL)},1,2);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);

                if isfield(problem,'jacobXdot')
                    Jxdot = feval(problem.jacobXdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:2
                                Jxdot{n} = sum(hbm.nonlin.hbm.Jxdot{n}.*States.df_dxdot(ijacobx,ijacobxnl,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Jxdot = repmat({zeros(NDofTot,NNLTot)},1,2);
                            for n = 1:2
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jxdot{n}(1:NDof,(2*l-1)*NNL+(1:NNL)) =  sum((-kHarm(l+1,n)*States.df_dxdot(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                                    Jxdot{n}(1:NDof,(2*l-0)*NNL+(1:NNL)) =  sum((-kHarm(l+1,n)*States.df_dxdot(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                    for l = 1:(NFreq-1)
                                        [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                        Jxdot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) =  2*sum((- kHarm(l+1,n)*States.df_dxdot(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxdot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) =  2*sum((- kHarm(l+1,n)*States.df_dxdot(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxdot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) = -2*sum((- kHarm(l+1,n)*States.df_dxdot(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                        Jxdot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) = -2*sum((- kHarm(l+1,n)*States.df_dxdot(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                                Jxdot{n} = Jxdot{n}(iRetain,iRetainNL);
                            end
                    end
                end
            end
            varargout{o} = Jxdot;
        case 'jacobUdot' %df_dU = w*Judot*U
            if ~hbm.dependence.udot
                Judot = repmat({zeros(NRetain,NInputTot)},1,2);
            else
                States.df_dudot = hbm_derivatives('nl','udot',States,hbm,problem);
                
                if isfield(problem,'jacobUdot')
                    Judot = feval(problem.jacobUdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:2
                                Judot{n} = sum(hbm.nonlin.hbm.Judot{n}.*States.df_dudot(ijacobx,ijacobu,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Judot = repmat({zeros(NDofTot,NInputTot)},1,2);
                            for n = 1:2
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Judot{n}(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum((-l*States.df_dudot.*sl)/(Nfft(1)*Nfft(2)),3);
                                    Judot{n}(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum((-l*States.df_dudot.*cl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                    for l = 1:(NFreq-1)
                                        [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                        Judot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,n)*States.df_dudot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Judot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,n)*States.df_dudot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Judot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,n)*States.df_dudot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                        Judot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,n)*States.df_dudot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                                Judot{n} = Judot{n}(iRetain,:);
                            end
                    end
                end
            end
            varargout{o} = Judot;
        case 'jacobXddot' %df_dX = w^2*Jxddot*X
            if ~hbm.dependence.xddot
                Jxddot = repmat({zeros(NRetain,NRetainNL)},1,3);
            else
                States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
                
                if isfield(problem,'jacobXdot')
                    Jxddot = feval(problem.jacobXdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:3
                                Jxddot{n} = sum(hbm.nonlin.hbm.Jxddot{n}.*States.df_dxddot(ijacobx,ijacobxnl,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Jxddot = repmat({zeros(NDof*NComp)},1,2);
                            for n = 1:2
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jxddot{n}(1:NDof,(2*l-1)*NNL+(1:NNL)) =  sum((-kHarm(l+1,n)^2*States.df_dxddot(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}(1:NDof,(2*l-0)*NNL+(1:NNL)) =  sum(( kHarm(l+1,n)^2*States.df_dxddot(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                    for l = 1:(NFreq-1)
                                        [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                        Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) =  2*sum((- kHarm(l+1,n)^2*States.df_dxddot(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) =  2*sum((  kHarm(l+1,n)^2*States.df_dxddot(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) = -2*sum((- kHarm(l+1,n)^2*States.df_dxddot(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                        Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) = -2*sum((  kHarm(l+1,n)^2*States.df_dxddot(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                                Jxddot{n} = Jxddot{n}(iRetain,iRetainNL);
                            end
                            
                            n = 3;
                            for l = 1:(NFreq-1)
                                [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                Jxddot{n}(1:NDof,(2*l-1)*NNL+(1:NNL)) =  sum((-prod(kHarm(l+1,:))*States.df_dxddot(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                                Jxddot{n}(1:NDof,(2*l-0)*NNL+(1:NNL)) =  sum(( prod(kHarm(l+1,:))*States.df_dxddot(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) =  2*sum((- prod(kHarm(l+1,:))*States.df_dxddot(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) =  2*sum((  prod(kHarm(l+1,:))*States.df_dxddot(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NNL+(1:NNL)) = -2*sum((- prod(kHarm(l+1,:))*States.df_dxddot(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NNL+(1:NNL)) = -2*sum((  prod(kHarm(l+1,:))*States.df_dxddot(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                            Jxddot{n} = Jxddot{n}(iRetain,iRetainNL);
                    end
                end
            end
            varargout{o} = Jxddot;
        case 'jacobUddot' %df_dU = w^2*Juddot*U
            if ~hbm.dependence.uddot
                Juddot = repmat({zeros(NRetain,NInputTot)},1,3);
            else
                States.df_duddot = hbm_derivatives('nl','uddot',States,hbm,problem);
                
                if isfield(problem,'jacobUddot')
                    Juddot = feval(problem.jacobUddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:3
                                Juddot{n} = sum(hbm.nonlin.hbm.Juddot{n}.*States.df_duddot(ijacobx,ijacobu,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Juddot = repmat({zeros(NDofTot,NInputTot)},1,2);
                            for n = 1:2
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Juddot{n}(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum((-kHarm(l+1,n)^2*States.df_duddot.*cl)/(Nfft(1)*Nfft(2)),3);
                                    Juddot{n}(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum(( kHarm(l+1,n)^2*States.df_duddot.*sl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                    for l = 1:(NFreq-1)
                                        [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                        Juddot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,n)^2*States.df_duddot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Juddot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum((  kHarm(l+1,n)^2*States.df_duddot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Juddot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,n)^2*States.df_duddot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                        Juddot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum((  kHarm(l+1,n)^2*States.df_duddot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                                Juddot{n} = Juddot{n}(iRetain,:);
                            end
                            n = 3;
                            for l = 1:(NFreq-1)
                                [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                Juddot{n}(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum((-prod(kHarm(l+1,:))*States.df_duddot.*cl)/(Nfft(1)*Nfft(2)),3);
                                Juddot{n}(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum(( prod(kHarm(l+1,:))*States.df_duddot.*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Juddot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum((- prod(kHarm(l+1,:))*States.df_duddot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Juddot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum((  prod(kHarm(l+1,:))*States.df_duddot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Juddot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum((- prod(kHarm(l+1,:))*States.df_duddot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    Juddot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum((  prod(kHarm(l+1,:))*States.df_duddot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                            Juddot{n} = Juddot{n}(iRetain,:);
                    end
                end
            end
            varargout{o} = Juddot;
        case 'floquet1xdot'
           if ~hbm.dependence.xdot
                D1 = zeros(NRetain,NRetainNL);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);
                
                if isfield(problem,'floquet1xdot')
                    D1 = feval(problem.floquet1xdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            D1 = sum(hbm.nonlin.hbm.Jx.*States.df_dxdot(ijacobx,ijacobxnl,:),3);
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1;theta2];
                            
                            D1 = zeros(NRetain,NRetainNL);
                            D1(1:NDof,1:NNL) = mean(States.df_dxdot(:,iNL,:),3);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                                D1(1:NDof,(2*l-1)*NDof+(1:NNL)) =  sum(( States.df_dxdot(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                                D1(1:NDof,(2*l-0)*NDof+(1:NNL)) =  sum((-States.df_dxdot(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(mtimesx(kHarm(k+1,:),theta));
                                D1((2*k-1)*NDof+(1:NDof),1:NNL) =  2*sum(States.df_dxdot(:,iNL,:).*ck/(Nfft(1)*Nfft(2)),3);
                                D1((2*k-0)*NDof+(1:NDof),1:NNL) = -2*sum(States.df_dxdot(:,iNL,:).*sk/(Nfft(1)*Nfft(2)),3);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NNL)) =  2*sum(( States.df_dxdot(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NNL)) =  2*sum((-States.df_dxdot(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NNL)) = -2*sum(( States.df_dxdot(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NNL)) = -2*sum((-States.df_dxdot(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                    end
                    D1 = D1(iRetain,iRetainNL);
                end
            end
            varargout{o} = D1;
        case 'floquet1xddot'
           if ~hbm.dependence.xddot
                D1dd = repmat({zeros(NRetain,NRetainNL)},1,2);
            else
                States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
                
                if isfield(problem,'floquet1xddot')
                    D1dd = feval(problem.floquet1xddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:2
                                D1dd{n} = sum(hbm.nonlin.hbm.Jxdot{n}.*States.df_dxddot(ijacobx,ijacobxnl,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1;theta2];
                            
                            for n = 1:2
                                D1dd{n} = zeros(NRetain,NRetainNL);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                                    D1dd{n}(1:NDof,(2*l-1)*NDof+(1:NNL)) =  sum((-kHarm(l+1,n)*States.df_dxddot(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                                    D1dd{n}(1:NDof,(2*l-0)*NDof+(1:NNL)) =  sum((-kHarm(l+1,n)*States.df_dxddot(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck, sk] = cosAndSin(mtimesx(kHarm(k+1,:),theta));
                                    for l = 1:(NFreq-1)
                                        [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                        D1dd{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NNL)) =  2*sum((-kHarm(l+1,n)*States.df_dxddot(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        D1dd{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NNL)) =  2*sum((-kHarm(l+1,n)*States.df_dxddot(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        D1dd{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NNL)) = -2*sum((-kHarm(l+1,n)*States.df_dxddot(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                        D1dd{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NNL)) = -2*sum((-kHarm(l+1,n)*States.df_dxddot(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                                D1dd{n} = D1dd{n}(iRetain,iRetainNL);
                            end
                    end
                end
            end
            varargout{o} = D1dd;
       case 'floquet2'
           if ~hbm.dependence.xddot
               D2 = zeros(NRetain);
           else
               States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
               
               if isfield(problem,'floquet2')
                   D2 = feval(problem.floquet2,States,hbm,problem);
               else
                   switch hbm.options.jacob_method
                       case 'mat'
                           D2 = sum(hbm.nonlin.hbm.Jx.*States.df_dxddot(ijacobx,ijacobxnl,:),3);
                       case 'sum'
                           theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                           theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                           [theta1,theta2] = ndgrid(theta1,theta2);
                           theta1 = permute(theta1(:)',[1 3 2]);
                           theta2 = permute(theta2(:)',[1 3 2]);
                           theta = [theta1;theta2];
                            
                           D2 = zeros(NRetain,NRetainNL);
                           D2(1:NDof,1:NNL) = mean(States.df_dxddot(:,iNL,:),3);
                           for l = 1:(NFreq-1)
                               [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                               D2(1:NDof,(2*l-1)*NDof+(1:NNL)) =  sum(( States.df_dxddot(:,iNL,:).*cl)/(Nfft(1)*Nfft(2)),3);
                               D2(1:NDof,(2*l-0)*NDof+(1:NNL)) =  sum((-States.df_dxddot(:,iNL,:).*sl)/(Nfft(1)*Nfft(2)),3);
                           end
                           for k = 1:(NFreq-1)
                               [ck, sk] = cosAndSin(mtimesx(kHarm(k+1,:),theta));
                               D2((2*k-1)*NDof+(1:NDof),1:NNL) =  2*sum(States.df_dxddot(:,iNL,:).*ck/(Nfft(1)*Nfft(2)),3);
                               D2((2*k-0)*NDof+(1:NDof),1:NNL) = -2*sum(States.df_dxddot(:,iNL,:).*sk/(Nfft(1)*Nfft(2)),3);
                               for l = 1:(NFreq-1)
                                   [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                                   D2((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NNL)) =  2*sum(( States.df_dxddot(:,iNL,:).*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                   D2((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NNL)) =  2*sum((-States.df_dxddot(:,iNL,:).*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                   D2((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NNL)) = -2*sum(( States.df_dxddot(:,iNL,:).*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                   D2((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NNL)) = -2*sum((-States.df_dxddot(:,iNL,:).*sl).*sk/(Nfft(1)*Nfft(2)),3);
                               end
                           end
                           
                   end
               end
            end
            varargout{o} = D2;
    end
end