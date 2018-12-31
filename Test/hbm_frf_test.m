function hbm_frf_test
problem = test_params;
NDof = problem.NDof;

hbm.harm.NHarm = 1;
hbm.harm.Nfft = 32;

hbm.cont.step0 = 0.1;
hbm.cont.max_step = 5;
hbm.cont.min_step = 1E-6;
hbm.cont.maxit = 100;

hbm.options.bAnalyticalDerivs = false;
hbm.options.bUseStandardHBM = true;
hbm.options.disp_dependence = true;
hbm.options.vel_dependence  = true;
hbm.options.freq_dependence = false;

hbm.options.cont_method = 'coco';
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

hbm = setuphbm(hbm,problem);

omega = sqrt(eig(problem.K,problem.M));
w0 = max(min(omega)-2,0.5);
wEnd = max(omega)+2;
A = 0.5;

tic;
hbm.pseudo.C = 1.08;
hbm.options.cont_method = 'pseudo';
[X,w,U,F,lambda] = hbm_frf(hbm,problem,w0,wEnd,A,[]);
X = permute((X(2,:,:)),[2 3 1]);
iStable = all(real(lambda)<0,1);
Xs = X; Xs(:,~iStable) = NaN;
Xus = X; Xus(:,iStable) = NaN;

t(1) = toc;
figure
hold on
col = lines(2);
for i = 1:2
    plot(w,abs(Xs(i,:)),'-','color',col(i,:));
    plot(w,abs(Xus(i,:)),'--','color',col(i,:));
end
return

tic;
hbm.arclength.C = 1;
hbm.options.cont_method = 'arclength';
[X2,w2] = hbm_frf(hbm,problem,w0,wEnd,A,[]);
X2 = permute((X2(2,:,:)),[2 3 1]);
t(2) = toc;

tic;
hbm.options.cont_method = 'coco';
[X3,w3] = hbm_frf(hbm,problem,w0,wEnd,A,[]);
X3 = permute((X3(2,:,:)),[2 3 1]);
t(3) = toc;

figure
ax(1) = subplot(2,1,1);
h = plot(w,abs(X),'r',w2,abs(X2),'g',w3,abs(X3),'b');
hold on
ylabel('|F| (mag)');

ax(2) = subplot(2,1,2);
plot(w,unwrap(angle(X),[],2),'r',w2,unwrap(angle(X2),[],2),'g',w3,unwrap(angle(X3),[],2),'b');
hold on
xlabel('\omega (rads)');
ylabel('\angle F (deg)');
linkaxes(ax,'x')
drawnow
legend(h(1:2:end),'pseudo','arclength','auto')

return

tic;
y = zeros(1,2*NDof);
w3 = [linspace(w0,wEnd,20)];% linspace(wEnd,w0,50)];
NCycle = 50;
Nfft = 101;
for i = 1:length(w3)
    y0 = y(end,:);
    U = A*test_model('excite',[],[],[],[],[],[],hbm,problem,w3(i));
    T = 2*pi/w3(i);
    [t,y] = ode45(@(t,y)test_odefun(t,y,w3(i),U,problem),[0 NCycle*T],y0);
    
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

