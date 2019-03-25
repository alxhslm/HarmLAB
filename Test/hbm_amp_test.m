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
hbm.cont.method = 'none';
sol1 = hbm_amp(hbm,problem,w0,A0,[],AEnd,[]);
X1 = cat(3,sol1.X);
X1 = permute(X1(2,:,:),[2 3 1]);
L = cat(2,sol1.L);
iStable = all(real(L)<0,1);
Xs = X1; Xs(:,~iStable) = NaN;
Xus = X1; Xus(:,iStable) = NaN;
A1 = cat(2,sol1.A);
t(1) = toc;
figure
hold on
col = lines(2);
for i = 1:2
    plot(A1,abs(Xs(i,:)),'-','color',col(i,:));
    plot(A1,abs(Xus(i,:)),'--','color',col(i,:));
end
leg{1} = 'none';

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'pseudo';
hbm.cont.predcorr.bMoorePenrose = 0;
sol2 = hbm_amp(hbm,problem,w0,A0,[],AEnd,[]);
X2 = cat(3,sol2.X);
X2 = permute(X2(2,:,:),[2 3 1]);
A2 = cat(2,sol2.A);
t(2) = toc;
leg{end+1} = 'pseudo';

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'arclength';
sol3 = hbm_amp(hbm,problem,w0,A0,[],AEnd,[]);
X3 = cat(3,sol3.X);
X3 = permute(X3(2,:,:),[2 3 1]);
A3 = cat(2,sol3.A);
t(3) = toc;
leg{end+1} = 'arclength';

figure
ax(1) = subplot(2,1,1);
h = plot(A1,abs(X1),'r',A2,abs(X2),'g',A3,abs(X3),'b');
hold on
ylabel('|F| (mag)');

ax(2) = subplot(2,1,2);
plot(A1,unwrap(angle(X1),[],2),'r',A2,unwrap(angle(X2),[],2),'g',A3,unwrap(angle(X3),[],2),'b');
hold on
xlabel('A (-)');
ylabel('\angle F (deg)');
linkaxes(ax,'x')
drawnow
legend(h(1:2:end),leg)

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