function varargout = hbm_nonlinear3d(command,hbm,problem,w0,Xp,Up)
if hbm.options.bUseStandardHBM
    [varargout{1:nargout}] = hbm_nonlinear(command,hbm,problem,w0(1),Xp,Up);
    return;
end

NInput = problem.NInput;
NDof = problem.NDof;

NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
NDofTot = NComp*NDof;
NInputTot = NComp*NInput;

%work out the time domain
X = unpackdof(Xp,NFreq-1,NDof);
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

Fp = packdof(F);

if ~iscell(command)
    command = {command};
end

varargout = cell(1,length(command));

ijacobx = hbm.nonlin.hbm.ijacobx;
ijacobu = hbm.nonlin.hbm.ijacobu;

for o = 1:length(command)
    switch command{o}
        case 'func'
            %all done
            varargout{o} = Fp;
        case 'derivW' %df_dw
            if ~hbm.dependence.w
                Dw = repmat({zeros(NDofTot,1)},1,2);
            else
                dfhbm_dw = hbm_derivatives('nl','w',States,hbm,problem);
                if isfield(problem,'derivW')
                    Dw = feval(problem.derivW,dfhbm_dw,hbm,problem);
                else
                    for n = 1:2
                        Dw{n} = packdof(hbm.nonlin.FFT*dfhbm_dw{n});
                    end
                end
            end
            varargout{o} = Dw;
        case 'jacobX'  %df_dX = Jx*X
            if ~hbm.dependence.x
                Jx = zeros(NDofTot,NDofTot);
            else
                States.df_dx = hbm_derivatives('nl','x',States,hbm,problem);
                
                if isfield(problem,'jacobX')
                    Jx = feval(problem.jacobX,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jx = sum(hbm.nonlin.hbm.Jx.*States.df_dx(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Jx = zeros(NDofTot,NDofTot);
                            Jx(1:NDof,1:NDof) = mean(States.df_dx,3);
                            for l = 1:(NFreq-1)
                                [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                Jx(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum( (States.df_dx.*cl)/(Nfft(1)*Nfft(2)),3);
                                Jx(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(-(States.df_dx.*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                Jx((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(States.df_dx.*ck/(Nfft(1)*Nfft(2)),3);
                                Jx((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(States.df_dx.*sk/(Nfft(1)*Nfft(2)),3);
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum( (States.df_dx.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum(-(States.df_dx.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum( (States.df_dx.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum(-(States.df_dx.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                    end
                end
            end
            varargout{o} = Jx;
        case 'jacobU' %df_dU = Ju*U
            if ~hbm.dependence.u
                Ju = zeros(NDofTot,NInputTot);
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
                    end
                end
            end
            varargout{o} = Ju;
        case 'jacobXdot' %df_dX = w*Jxdot*X
            if ~hbm.dependence.xdot
                Jxdot = repmat({zeros(NDofTot,NDofTot)},1,2);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);

                if isfield(problem,'jacobXdot')
                    Jxdot = feval(problem.jacobXdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:2
                                Jxdot{n} = sum(hbm.nonlin.hbm.Jxdot{n}.*States.df_dxdot(ijacobx,ijacobx,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1; theta2];
                            
                            Jxdot = repmat({zeros(NDofTot,NDofTot)},1,2);
                            for n = 1:2
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jxdot{n}(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,n)*States.df_dxdot.*sl)/(Nfft(1)*Nfft(2)),3);
                                    Jxdot{n}(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,n)*States.df_dxdot.*cl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                    for l = 1:(NFreq-1)
                                        [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                        Jxdot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,n)*States.df_dxdot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxdot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,n)*States.df_dxdot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxdot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,n)*States.df_dxdot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                        Jxdot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,n)*States.df_dxdot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                            end
                    end
                end
            end
            varargout{o} = Jxdot;
        case 'jacobUdot' %df_dU = w*Judot*U
            if ~hbm.dependence.udot
                Judot = repmat({zeros(NDofTot,NInputTot)},1,2);
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
                            end
                    end
                end
            end
            varargout{o} = Judot;
        case 'jacobXddot' %df_dX = w^2*Jxddot*X
            if ~hbm.dependence.xddot
                Jxddot = repmat({zeros(NDofTot,NDofTot)},1,3);
            else
                States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
                
                if isfield(problem,'jacobXdot')
                    Jxddot = feval(problem.jacobXdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:3
                                Jxddot{n} = sum(hbm.nonlin.hbm.Jxddot{n}.*States.df_dxddot(ijacobx,ijacobx,:),3);
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
                                    Jxddot{n}(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,n)^2*States.df_dxddot.*cl)/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(( kHarm(l+1,n)^2*States.df_dxddot.*sl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                    for l = 1:(NFreq-1)
                                        [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                        Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,n)^2*States.df_dxddot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((  kHarm(l+1,n)^2*States.df_dxddot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,n)^2*States.df_dxddot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                        Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((  kHarm(l+1,n)^2*States.df_dxddot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                            end
                            
                            n = 3;
                            for l = 1:(NFreq-1)
                                [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                Jxddot{n}(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-prod(kHarm(l+1,:))*States.df_dxddot.*cl)/(Nfft(1)*Nfft(2)),3);
                                Jxddot{n}(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(( prod(kHarm(l+1,:))*States.df_dxddot.*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck,sk] = cosAndSin(kHarm(k+1,1)*theta(1,:,:) + kHarm(k+1,2)*theta(2,:,:));
                                for l = 1:(NFreq-1)
                                    [cl,sl] = cosAndSin(kHarm(l+1,1)*theta(1,:,:) + kHarm(l+1,2)*theta(2,:,:));
                                    Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- prod(kHarm(l+1,:))*States.df_dxddot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((  prod(kHarm(l+1,:))*States.df_dxddot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- prod(kHarm(l+1,:))*States.df_dxddot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    Jxddot{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((  prod(kHarm(l+1,:))*States.df_dxddot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                    end
                end
            end
            varargout{o} = Jxddot;
        case 'jacobUddot' %df_dU = w^2*Juddot*U
            if ~hbm.dependence.uddot
                Juddot = repmat({zeros(NDofTot,NInputTot)},1,3);
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
                    end
                end
            end
            varargout{o} = Juddot;
        case 'floquet1xdot'
           if ~hbm.dependence.xdot
                D1 = zeros(NDofTot,NDofTot);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);
                
                if isfield(problem,'floquet1xdot')
                    D1 = feval(problem.floquet1xdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            D1 = sum(hbm.nonlin.hbm.Jx.*States.df_dxdot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1;theta2];
                            
                            D1 = zeros(NDofTot,NDofTot);
                            D1(1:NDof,1:NDof) = mean(States.df_dxdot,3);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                                D1(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum(( States.df_dxdot.*cl)/(Nfft(1)*Nfft(2)),3);
                                D1(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-States.df_dxdot.*sl)/(Nfft(1)*Nfft(2)),3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(mtimesx(kHarm(k+1,:),theta));
                                D1((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(States.df_dxdot.*ck/(Nfft(1)*Nfft(2)),3);
                                D1((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(States.df_dxdot.*sk/(Nfft(1)*Nfft(2)),3);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum(( States.df_dxdot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-States.df_dxdot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum(( States.df_dxdot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-States.df_dxdot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                end
                            end
                    end
                end
            end
            varargout{o} = D1;
        case 'floquet1xddot'
           if ~hbm.dependence.xddot
                D1dd = repmat({zeros(NDofTot,NDofTot)},1,2);
            else
                States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
                
                if isfield(problem,'floquet1xddot')
                    D1dd = feval(problem.floquet1xddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            for n = 1:2
                                D1dd{n} = sum(hbm.nonlin.hbm.Jxdot{n}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                            end
                        case 'sum'
                            theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                            theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                            [theta1,theta2] = ndgrid(theta1,theta2);
                            theta1 = permute(theta1(:)',[1 3 2]);
                            theta2 = permute(theta2(:)',[1 3 2]);
                            theta = [theta1;theta2];
                            
                            for n = 1:2
                                D1dd{n} = zeros(NDofTot,NDofTot);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                                    D1dd{n}(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,n)*States.df_dxddot.*sl)/(Nfft(1)*Nfft(2)),3);
                                    D1dd{n}(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,n)*States.df_dxddot.*cl)/(Nfft(1)*Nfft(2)),3);
                                end
                                for k = 1:(NFreq-1)
                                    [ck, sk] = cosAndSin(mtimesx(kHarm(k+1,:),theta));
                                    for l = 1:(NFreq-1)
                                        [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                        D1dd{n}((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,n)*States.df_dxddot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                        D1dd{n}((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,n)*States.df_dxddot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                        D1dd{n}((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,n)*States.df_dxddot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                                        D1dd{n}((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,n)*States.df_dxddot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                    end
                                end
                            end
                    end
                end
            end
            varargout{o} = D1dd;
       case 'floquet2'
           if ~hbm.dependence.xddot
               D2 = zeros(NDofTot);
           else
               States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
               
               if isfield(problem,'floquet2')
                   D2 = feval(problem.floquet2,States,hbm,problem);
               else
                   switch hbm.options.jacob_method
                       case 'mat'
                           D2 = sum(hbm.nonlin.hbm.Jx.*States.df_dxddot(ijacobx,ijacobx,:),3);
                       case 'sum'
                           theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
                           theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
                           [theta1,theta2] = ndgrid(theta1,theta2);
                           theta1 = permute(theta1(:)',[1 3 2]);
                           theta2 = permute(theta2(:)',[1 3 2]);
                           theta = [theta1;theta2];
                            
                           D2 = zeros(NDofTot,NDofTot);
                           D2(1:NDof,1:NDof) = mean(States.df_dxddot,3);
                           for l = 1:(NFreq-1)
                               [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                               D2(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum(( States.df_dxddot.*cl)/(Nfft(1)*Nfft(2)),3);
                               D2(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-States.df_dxddot.*sl)/(Nfft(1)*Nfft(2)),3);
                           end
                           for k = 1:(NFreq-1)
                               [ck, sk] = cosAndSin(mtimesx(kHarm(k+1,:),theta));
                               D2((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(States.df_dxddot.*ck/(Nfft(1)*Nfft(2)),3);
                               D2((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(States.df_dxddot.*sk/(Nfft(1)*Nfft(2)),3);
                               for l = 1:(NFreq-1)
                                   [cl, sl] = cosAndSin(mtimesx(kHarm(l+1,:),theta));
                                   D2((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum(( States.df_dxddot.*cl).*ck/(Nfft(1)*Nfft(2)),3);
                                   D2((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-States.df_dxddot.*sl).*ck/(Nfft(1)*Nfft(2)),3);
                                   D2((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum(( States.df_dxddot.*cl).*sk/(Nfft(1)*Nfft(2)),3);
                                   D2((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-States.df_dxddot.*sl).*sk/(Nfft(1)*Nfft(2)),3);
                               end
                           end
                           
                   end
               end
            end
            varargout{o} = D2;
    end
end