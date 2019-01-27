function hbm_amp_test
problem = test_params;
NDof = problem.NDof;

hbm.harm.NHarm = 2;
hbm.harm.Nfft = 32;

hbm.options.bUseStandardHBM = true;
hbm.dependence.x = true;
hbm.dependence.xdot  = true;
hbm.dependence.w = false;

hbm.scaling.tol = 1;

hbm.cont.step0 = 1E-2;
hbm.cont.max_step = 1E-1;
hbm.cont.min_step = 1E-6;

hbm.cont.method = 'predcorr';
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

[hbm,problem] = setuphbm(hbm,problem);

omega = sqrt(eig(problem.K,problem.M));
w0 = omega(1);
wEnd = max(omega)+2;
A0 = 10;
AEnd = 30;

tic;
hbm.cont.method = 'coco';
sol1 = hbm_amp(hbm,problem,w0,A0,[],AEnd,[]);
X = permute((sol1.X(2,:,:)),[2 3 1]);
iStable = all(real(sol1.L)<0,1);
Xs = X; Xs(:,~iStable) = NaN;
Xus = X; Xus(:,iStable) = NaN;

t(1) = toc;
figure
hold on
col = lines(2);
for i = 1:2
    plot(sol1.A,abs(Xs(i,:)),'-','color',col(i,:));
    plot(sol1.A,abs(Xus(i,:)),'--','color',col(i,:));
end

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'pseudo';
hbm.cont.predcorr.bMoorePenrose = 0;
sol2 = hbm_amp(hbm,problem,w0,A0,[],AEnd,[]);
X2 = permute((sol2.X(2,:,:)),[2 3 1]);
t(2) = toc;

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'arclength';
sol3 = hbm_amp(hbm,problem,w0,A0,[],AEnd,[]);
X3 = permute((sol3.X(2,:,:)),[2 3 1]);
t(3) = toc;

figure
ax(1) = subplot(2,1,1);
h = plot(sol1.A,abs(X),'r',sol2.A,abs(X2),'g',sol3.A,abs(X3),'b');
hold on
ylabel('|F| (mag)');

ax(2) = subplot(2,1,2);
plot(sol1.A,unwrap(angle(X),[],2),'r',sol2.A,unwrap(angle(X2),[],2),'g',sol3.A,unwrap(angle(X3),[],2),'b');
hold on
xlabel('A (-)');
ylabel('\angle F (deg)');
linkaxes(ax,'x')
drawnow
legend(h(1:2:end),'pseudo','arclength','coco')

return

tic;
y = zeros(1,2*NDof);
sol3.w = [linspace(w0,wEnd,20)];% linspace(wEnd,w0,50)];
NCycle = 50;
Nfft = 101;
for i = 1:length(sol3.w)
    y0 = y(end,:);
    U = A*test_model('excite',[],[],[],[],[],[],hbm,problem,sol3.w(i));
    T = 2*pi/sol3.w(i);
    [t,y] = ode45(@(t,y)test_odefun(t,y,sol3.w(i),U,problem),[0 NCycle*T],y0);
    
    ii = find(t>(t(end)-T),1);
    t = t(ii:end) - t(ii);
    y = y(ii:end,:);
    
    tReg = linspace(0,T,Nfft)';
    Fs = Nfft/T;
    yReg = interp1(t,y,tReg,'linear','extrap');
    
    Yfft = 2*fft(yReg,[],1)/Nfft;
    Wfft = (0:(Nfft-1))'/Nfft * 2*pi* Fs;
        
    Xode(i,:) = Yfft(2,1:NDof);
    wode(i) = Wfft(2);
end
toc;

plot(ax(1),wode,abs(Xode));
plot(ax(2),wode,unwrap(angle(Xode),[],2));