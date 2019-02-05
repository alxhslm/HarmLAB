function hbm_frf_plot(command,hbm,problem,results)
persistent fig hSuccess hWarn hErr X W
if hbm.cont.bUpdate
    switch command
        case 'init'
            if ~isempty(fig) && ishandle(fig)
                close(fig)
            end
                          
            if isempty(results)
                Xi = zeros(hbm.harm.NFreq,problem.NDof);
                wi = NaN;
            else
                Xi = results.X;
                wi = results.w;
            end
            W = getfrequencies(wi,hbm);
            X = Xi;
            
            [xlin, wlin] = getLinearReponse(hbm,problem,X,wi,results.A);
            [fig,hSuccess,hWarn,hErr] = createFRF(hbm,problem,X,W,xlin,wlin);

        case {'data','err','warn'}
            if ~ishandle(fig(1))
                [xlin, wlin] = getLinearReponse(hbm,problem,X(:,:,1),W(1),results.A);
                [fig,hSuccess,hWarn,hErr] = createFRF(hbm,problem,X,W,xlin,wlin);
            end
            
            X(:,:,end+1) = results.X;
            W(:,end+1) = getfrequencies(results.w,hbm);
            Xabs = abs(X); Xabs = Xabs(:,:,end);
            Xph = unwrap(angle(X)); Xph = Xph(:,:,end);
            Wfreq = W(:,end);
            if any(strcmpi(command,{'data','warn'}))
                update_handles(hSuccess,Xabs,Xph,Wfreq,hbm,problem)
                %update our progress
                
                if strcmpi(command,'warn')
                    %warning, overlay in blue
                    update_handles(hWarn,Xabs,Xph,Wfreq,hbm,problem)
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
                W(:,end) = [];
                update_handles(hErr,Xabs,Xph,Wfreq,hbm,problem)
            end

            drawnow
        case 'close'
            close(fig)
            hSuccess = [];
            hWarn = [];
            hErr = [];
    end
end

function update_handles(han,Xabs,Xph,W,hbm,problem)
for i = 1:length(hbm.harm.iHarmPlot)
    for j = 1:length(problem.iDofPlot)
        w = [get(han{1}(i,j),'xdata'),W(hbm.harm.iHarmPlot(i),:)];
        mag = [get(han{1}(i,j),'ydata'),permute(Xabs(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:),[1 3 2])];
        ph  = [get(han{2}(i,j) ,'ydata'),permute(Xph(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:),[1 3 2])];
        set(han{1}(i,j),'xdata',w,'ydata',mag);
        set(han{2}(i,j) ,'xdata',w,'ydata',ph);
    end
end

function ws = getfrequencies(w0,hbm)
w = (hbm.harm.rFreqBase.*hbm.harm.rFreqRatio)'*w0;
ws = abs(hbm.harm.kHarm(:,1)*w(1,:) + hbm.harm.kHarm(:,2)*w(2,:));
ws(ws == 0) = w0;

function [fMag,hSuccess,hWarn,hErr] = createFRF(hbm,problem,x,w,xlin,wlin0)
matlabPos = getMatlabSize;
figPos = matlabPos;
figPos(4) = matlabPos(4)/2;
figPos(2) = matlabPos(2) + figPos(4);
fMag = figure('Name',[problem.name],'OuterPosition',figPos,'WindowStyle', 'Docked');

wlin = getfrequencies(wlin0,hbm);
wlim = getfrequencies([problem.wMin problem.wMax],hbm);

for i = 1:length(hbm.harm.iHarmPlot)
    for j = 1:length(problem.iDofPlot)
        tmp = subplot(length(problem.iDofPlot),length(hbm.harm.iHarmPlot),(j-1)*length(hbm.harm.iHarmPlot) + i,'Parent',fMag);        
        Xij = squeeze(xlin(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:));
        wij = wlin(hbm.harm.iHarmPlot(i),:);
        [tmp2,hLin{1}(i,j),hLin{2}(i,j)] = plotyy(tmp,wij,abs(Xij),wij,unwrap(angle(Xij)));
        ax{1}(i,j) = tmp2(1); ax{2}(i,j) = tmp2(2);
        hold(tmp2(1), 'on');
        hold(tmp2(2), 'on');
    end
end

for i = 1:length(hbm.harm.iHarmPlot)
    for j = 1:length(problem.iDofPlot)
        hSuccess{1}(i,j)  = plot(ax{1}(i,j),w(hbm.harm.iHarmPlot(i),:),abs(squeeze(x(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:))),'g.-');
        hSuccess{2}(i,j)  = plot(ax{2}(i,j),w(hbm.harm.iHarmPlot(i),:),unwrap(angle(squeeze(x(hbm.harm.iHarmPlot(i),problem.iDofPlot(j),:)))),'g.-');
        
        for k = 1:2
            hWarn{k}(i,j) = plot(ax{k}(i,j),NaN,NaN,'b.');
            hErr{k}(i,j)  = plot(ax{k}(i,j),NaN,NaN,'r.');
        
            %xlim(ax{k}(i,j),wlim(hbm.harm.iHarmPlot(i),:));
            set(ax{k}(i,j),'XLimMode','auto')
            set(ax{k}(i,j),'YLimMode','auto')
        end
        
        if j==length(problem.iDofPlot)
            xlabel(ax{1}(i,j),'\omega (rads)')
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

function [xlin, wlin] = getLinearReponse(hbm,problem,X,w0,A)
%find the linearised contribution to the stiffness/damping due from the non-linearity
w0 = w0*hbm.harm.rFreqRatio;
wB = w0.*hbm.harm.rFreqBase;
NHarm = hbm.harm.NHarm;
Nfft = hbm.harm.Nfft;

w = hbm.harm.kHarm*wB';
x0 = X(1,:).';
U = A*feval(problem.excite,hbm,problem,w0);
u0 = U(1,:).';

States = hbm_states3d(w0,X,U,hbm);

States.f = feval(problem.model,'nl',States,hbm,problem);

[K_nl, C_nl, M_nl]  = hbm_derivatives('nl',{'x','xdot','xddot'},States,hbm,problem);
[Ku_nl,Cu_nl,Mu_nl] = hbm_derivatives('nl',{'u','udot','uddot'},States,hbm,problem);

K_nl  = mean(K_nl,3);  C_nl  = mean(C_nl,3);  M_nl  = mean(M_nl,3); 
Ku_nl = mean(Ku_nl,3); Cu_nl = mean(Cu_nl,3); Mu_nl = mean(Mu_nl,3);

M  = problem.M + M_nl;
C  = problem.C + C_nl;
K  = problem.K + K_nl;

Mu = problem.Mu + Mu_nl;
Cu = problem.Cu + Cu_nl;
Ku = problem.Ku + Ku_nl;

%now loop over all the frequencies
wlin = linspace(problem.wMin,problem.wMax,1000);
xlin = zeros(hbm.harm.NFreq,problem.NDof,length(wlin));

NFreq = hbm.harm.NFreq;

%work out frequencies
w = getfrequencies(wlin,hbm);
for i = 1:length(wlin)
    w0 = wlin(i);
    
    U = A*feval(problem.excite,hbm,problem,w0);
    
    for k = 1:NFreq
        Fe = (Ku + 1i*w(k,i)*Cu - w(k,i)^2*Mu)*U(k,:).';
        H = K + 1i*w(k,i)*C - w(k,i)^2 * M;
        xlin(k,:,i) = (H\Fe).';
    end
end
xlin(1,:,:) = xlin(1,:,:) + x0.';