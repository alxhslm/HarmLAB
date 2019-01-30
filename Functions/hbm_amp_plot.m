function hbm_amp_plot(command,hbm,problem,x,a)
persistent fig hSuccess hWarn hErr X A
w0 = problem.w0;

if hbm.cont.bUpdate
    switch command
        case 'init'
            if ~isempty(fig) && ishandle(fig)
                close(fig)
            end
                          
            if isempty(x) || isempty(A)
                X = zeros(hbm.harm.NFreq,problem.NDof);
                A = NaN;
            else
                X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
                A = a;
            end
            [xlin, Alin] = getLinearReponse(hbm,problem,X,w0);
            [fig,hSuccess,hWarn,hErr] = createFRF(hbm,problem,X,A,xlin,Alin);

        case {'data','err','warn'}
            if ~ishandle(fig(1))
                [xlin, Alin] = getLinearReponse(hbm,problem,X(:,:,1),w0);
                [fig,hSuccess,hWarn,hErr] = createFRF(hbm,problem,X,A,xlin,Alin);
            end
            
            X(:,:,end+1) = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
            A(:,end+1) = a;
            Xabs = abs(X); Xabs = Xabs(:,:,end);
            Xph = unwrap(angle(X)); Xph = Xph(:,:,end);
            if any(strcmpi(command,{'data','warn'}))
                update_handles(hSuccess,Xabs,Xph,a,hbm,problem)
                %update our progress
                
                if strcmpi(command,'warn')
                    %warning, overlay in blue
                    update_handles(hWarn,Xabs,Xph,a,hbm,problem)
                end
                
                %reset the error points
                for i = 1:length(hbm.harm.iHarmPlot)
                    for j = 1:length(problem.iDofPlot)
                        for k = 1:2
                            set(hWarn{k}(i,j),'xdata',NaN,'ydata',NaN);
                            set(hErr{k}(i,j) ,'xdata',NaN,'ydata',NaN);
                        end
                    end
                end
            else
                %error
                X(:,:,end) = [];
                A(:,end) = [];
                update_handles(hErr,Xabs,Xph,a,hbm,problem)
            end

            drawnow
        case 'close'
            close(fig)
            hSuccess = [];
            hWarn = [];
            hErr = [];
    end
end

function update_handles(han,Xabs,Xph,A,hbm,problem)
for i = 1:length(hbm.harm.iHarmPlot)
    for j = 1:length(problem.iDofPlot)
        a = [get(han{1}(i,j),'xdata'),A];
        mag = [get(han{1}(i,j),'ydata'),permute(Xabs(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:),[1 3 2])];
        ph  = [get(han{2}(i,j) ,'ydata'),permute(Xph(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:),[1 3 2])];
        set(han{1}(i,j),'xdata',a,'ydata',mag);
        set(han{2}(i,j) ,'xdata',a,'ydata',ph);
    end
end

function [fMag,hSuccess,hWarn,hErr] = createFRF(hbm,problem,x,A,xlin,Alin)
% matlabPos = getMatlabSize;
% figPos = matlabPos;
% figPos(4) = matlabPos(4)/2;
% figPos(2) = matlabPos(2) + figPos(4);
fMag = figure('Name',[problem.name]);%,'OuterPosition',figPos,'WindowStyle', 'Docked');

for i = 1:length(hbm.harm.iHarmPlot)
    for j = 1:length(problem.iDofPlot)
        tmp = subplot(length(problem.iDofPlot),length(hbm.harm.iHarmPlot),(j-1)*length(hbm.harm.iHarmPlot) + i,'Parent',fMag);        
        Xij = squeeze(xlin(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:));
        [tmp2,hLin{1}(i,j),hLin{2}(i,j)] = plotyy(tmp,Alin,abs(Xij),Alin,unwrap(angle(Xij)));
        ax{1}(i,j) = tmp2(1); ax{2}(i,j) = tmp2(2);
        hold(tmp2(1), 'on');
        hold(tmp2(2), 'on');
    end
end

for i = 1:length(hbm.harm.iHarmPlot)
    for j = 1:length(problem.iDofPlot)
        hSuccess{1}(i,j)  = plot(ax{1}(i,j),A,abs(squeeze(x(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:))),'g.-');
        hSuccess{2}(i,j)  = plot(ax{2}(i,j),A,unwrap(angle(squeeze(x(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:)))),'g.-');
        
        for k = 1:2
            hWarn{k}(i,j) = plot(ax{k}(i,j),NaN,NaN,'b.');
            hErr{k}(i,j)  = plot(ax{k}(i,j),NaN,NaN,'r.');
        
            %xlim(ax{k}(i,j),wlim(hbm.harm.iHarmPlot(i),:));
            set(ax{k}(i,j),'XLimMode','auto')
            set(ax{k}(i,j),'YLimMode','auto')
        end
        
        if j==length(problem.iDofPlot)
            xlabel(ax{1}(i,j),'A (-)')
        end
        if j==1
            title(ax{1}(i,j),harmonicName(hbm,i))
        end
        
        if i == 1
            ylabel(ax{1}(i,j),sprintf('|Dof #%d|',problem.iDofPlot(j)))
        end
        
        if i == length(hbm.harm.iHarmPlot)
            ylabel(ax{2}(i,j),sprintf('\\angle Dof #%d',problem.iDofPlot(j)))
        end
    end
end

function s = harmonicName(hbm,i)
k1 = hbm.harm.kHarm(hbm.harm.iHarmPlot(i),1)*hbm.harm.rFreqBase(1);
k2 = hbm.harm.kHarm(hbm.harm.iHarmPlot(i),2)*hbm.harm.rFreqBase(2);
if k1 == 0
    s1 = '';
elseif k1 == 1
    s1 = '\omega_1';
else
    s1 = sprintf('%s\\omega_1',num2frac(k1));
end
if k2 == 0
    s2 = '';
elseif k2 == 1
    if ~isempty(s1)
        s2 = '+\omega_2';
    else
        s2 = '\omega_2';
    end
elseif k2 == -1
        s2 = '-\omega_2';
else
    if ~isempty(s1)
        s2 = sprintf('%+s\\omega_2',num2frac(k2));
    else
        s2 = sprintf('%s\\omega_2',num2frac(k2));
    end
end

if isempty(s2)
    s = s1;
elseif isempty(s1)
    s = s2;
else
    s = [s1 ' ' s2];
end

if isempty(s)
    s = '0';
end

function [xlin, Alin] = getLinearReponse(hbm,problem,X,w0)
%find the linearised contribution to the stiffness/damping due from the non-linearity
w0 = w0*hbm.harm.rFreqRatio;
wB = w0.*hbm.harm.rFreqBase;
NHarm = hbm.harm.NHarm;
Nfft = hbm.harm.Nfft;

Alin = linspace(problem.A0,problem.AEnd,1000);

w = hbm.harm.kHarm*wB';
x0 = X(1,:).';
U = feval(problem.excite,hbm,problem,w0);

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%precompute the external inputs
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

x     = freq2time3d(X,NHarm,hbm.harm.iSub,Nfft);
xdot  = freq2time3d(Xdot,NHarm,hbm.harm.iSub,Nfft);
xddot = freq2time3d(Xddot,NHarm,hbm.harm.iSub,Nfft);

%create the vector of inputs
u     = freq2time3d(U,NHarm,hbm.harm.iSub,Nfft);
udot  = freq2time3d(Udot,NHarm,hbm.harm.iSub,Nfft);
uddot = freq2time3d(Uddot,NHarm,hbm.harm.iSub,Nfft);

% %create the time series from the fourier series
% x     = repmat(x0,1,prod(Nfft))';
% xdot  = 0*x;
% xddot = 0*x;
% 
% %create the vector of inputs
% u     = repmat(u0,1,prod(Nfft))';
% udot  = 0*u;
% uddot = 0*u;

%work out the time vector
t1 = (0:Nfft(1)-1)/Nfft(1)*2*pi/wB(1);
t2 = (0:Nfft(2)-1)/Nfft(2)*2*pi/wB(2);
[t1,t2] = ndgrid(t1,t2);
t = [t1(:) t2(:)];

xlin = zeros(hbm.harm.NFreq,problem.NDof,length(Alin));

%now loop over all the amplitudes

for i = 1:length(Alin)
    f0 = feval(problem.model,'nl',t',x',xdot',xddot',Alin(i)*u',Alin(i)*udot',Alin(i)*uddot',hbm,problem,w0).';

    [K_nl, C_nl, M_nl]  = hbm_derivatives('nl',{'x','xdot','xddot'},t,x,xdot,xddot,Alin(i)*u,Alin(i)*udot,Alin(i)*uddot,f0,hbm,problem,w0);
    [Ku_nl,Cu_nl,Mu_nl] = hbm_derivatives('nl',{'u','udot','uddot'},t,x,xdot,xddot,Alin(i)*u,Alin(i)*udot,Alin(i)*uddot,f0,hbm,problem,w0);

    K_nl  = mean(K_nl,3);  C_nl  = mean(C_nl,3);  M_nl  = mean(M_nl,3); 
    Ku_nl = mean(Ku_nl,3); Cu_nl = mean(Cu_nl,3); Mu_nl = mean(Mu_nl,3);

    M  = problem.M + M_nl;
    C  = problem.C + C_nl;
    K  = problem.K + K_nl;

    Mu = problem.Mu + Mu_nl;
    Cu = problem.Cu + Cu_nl;
    Ku = problem.Ku + Ku_nl;

    NFreq = hbm.harm.NFreq;
    
    for k = 1:NFreq
        Fe = (Ku + 1i*w(k)*Cu - w(k)^2*Mu)*Alin(i)*U(k,:).';
        H = K + 1i*w(k)*C - w(k)^2 * M;
        xlin(k,:,i) = (H\Fe).';
    end
end

xlin(1,:,:) = xlin(1,:,:) + x0.';