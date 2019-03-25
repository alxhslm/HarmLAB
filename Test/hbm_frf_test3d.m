function hbm_frf_test3d
problem = test_params;
NDof = problem.NDof;

hbm.harm.rFreqRatio = 1;
hbm.harm.NHarm = 2;
hbm.harm.Nfft = 32;

if 1
    hbm.harm.rFreqRatio(end+1) = 1.376;
    hbm.harm.NHarm(end+1) = 1;
    hbm.harm.Nfft(end+1) = 16;
end

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

hbm.harm.iHarmPlot = [1 2 3];

[hbm,problem] = setuphbm(hbm,problem);

omega = sqrt(eig(problem.K,problem.M));
w0 = max(min(omega)-2,0.5);
wEnd = max(omega)+2;
A = 10;

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'pseudo';
sol1 = hbm_frf(hbm,problem,A,w0,[],wEnd,[]);
X1 = permute((sol1.X(2,:,:)),[2 3 1]);
iStable = all(real(sol1.L)<0,1);
Xs = X1; Xs(:,~iStable) = NaN;
Xus = X1; Xus(:,iStable) = NaN;
w1 = sol1.w;
t(1) = toc;
figure
hold on
col = lines(2);
for i = 1:2
    plot(w1,abs(Xs(i,:)),'-','color',col(i,:));
    plot(w1,abs(Xus(i,:)),'--','color',col(i,:));
end

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'arclength';
sol2 = hbm_frf(hbm,problem,A,w0,[],wEnd,[]);
X2 = permute((sol2.X(2,:,:)),[2 3 1]);
w2 = sol2.w;
t(2) = toc;

tic;
hbm.cont.method = 'coco';
sol3 = sol2;%hbm_frf(hbm,problem,A,w0,[],wEnd,[]);
X3 = permute((sol3.X(2,:,:)),[2 3 1]);
w3 = sol3.w;
t(3) = toc;

figure
ax(1) = subplot(2,1,1);
h = plot(w1,abs(X1),'r',w2,abs(X2),'g',w3,abs(X3),'b');
hold on
ylabel('|F| (mag)');

ax(2) = subplot(2,1,2);
plot(w1,unwrap(angle(X1),[],2),'r',w2,unwrap(angle(X2),[],2),'g',w3,unwrap(angle(X3),[],2),'b');
hold on
xlabel('\omega (rads)');
ylabel('\angle F (deg)');
linkaxes(ax,'x')
drawnow
legend(h(1:2:end),'pseudo','arclength','coco')

save(fullfile(onedriveroot,'MATLAB','master'),'w1','X1','w2','X2','w3','X3')

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