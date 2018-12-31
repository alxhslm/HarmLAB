function hbm_bb_test3d
problem = test_params;
P = problem.P;

hbm.harm.rFreqRatio = 1;
hbm.harm.NHarm = 1;
hbm.harm.Nfft = 16;
if 1
    hbm.harm.rFreqRatio(end+1) = 1.376;
    hbm.harm.NHarm(end+1) = 2;
    hbm.harm.Nfft(end+1) = 16;
        hbm.harm.kHarm = [0 0;
                          1 0;
                          0 1;
                          0 2];
end

hbm.cont.step0  = 1;
hbm.cont.min_step = 0.1;
hbm.cont.max_step = 100;
hbm.cont.maxit = 100;

hbm.options.disp_dependence = true;
hbm.options.vel_dependence  = true;
hbm.options.freq_dependence = false;

hbm.options.cont_method = 'pseudo';
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

hbm.options.bAnalyticalDerivs = false;
hbm = setuphbm(hbm,problem);

iMode = 1;
omega = sqrt(eig(problem.K,problem.M));
w0 = max(omega(iMode)-2,0.5);
wEnd = omega(iMode)+1.5;

figure
ax(1) = subplot(2,1,1);
ylabel('|F| (mag)');
hold on
ax(2) = subplot(2,1,2);
hold on
xlabel('\omega (rads)');
ylabel('\angle F (deg)');
linkaxes(ax,'x')
xlim([w0,wEnd]);

for A = [1 2 3]
    file = ['FRF_' sprintf('A = %0.1f',A) '.mat'];
    if ~isfile(file)
        [x,w,u] = hbm_frf(hbm,problem,w0,wEnd,A,[]);
        save(file,'x','w','u')
    else
        load(file)
    end
    X = permute((x(2,:,:)),[2 3 1]);
    U = permute((u(2,:,:)),[2 3 1]);
    H = X(P.iDof,:)./U(P.iInput,:);
    plot(ax(1),w,abs(H));
    plot(ax(2),w,unwrap(angle(H),[],2));
    
    [~,ii] = max(H);
    [Xres,wres,Ures] = hbm_resonance(hbm,problem,w(ii),A,x(:,:,ii));
    Xres = permute((Xres(2,:,:)),[2 3 1]);
    Ures = permute((Ures(2,:,:)),[2 3 1]);
    Hres = Xres(P.iDof,:)./Ures(P.iInput,:);
    plot(ax(1),wres,abs(Hres),'o')
    plot(ax(2),wres,angle(Hres),'o')
    drawnow
end

bb = hbm;
bb.options.cont_method = 'none';

tic;
A0 = 1;
Aend = 3;
vomit
[X_bb,w_bb,~,U_bb] = hbm_bb(bb,problem,A0,Aend,omega(iMode),[]);
X_bb = permute(X_bb(2,:,:),[2 3 1]);
U_bb = permute(U_bb(2,:,:),[2 3 1]);
H_bb = X_bb(P.iDof,:)./U_bb(P.iInput,:);
toc;
plot(ax(1),w_bb,abs(H_bb));
plot(ax(2),w_bb,unwrap(angle(H_bb),[],2));