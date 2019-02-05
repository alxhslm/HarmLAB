function hbm_bb_test
problem = test_params;

P = problem.P;

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

hbm.options.bAnalyticalDerivs = true;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'pseudo';
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

problem.res.iDof = P.iDof;
problem.res.iInput = P.iInput;
problem.res.iHarm = 2;
problem.res.sign = 1;
problem.res.output = 'x';
problem.res.input = 'u';

[hbm,problem] = setuphbm(hbm,problem);

iMode = 1;
omega = sqrt(eig(problem.K,problem.M));
w0 = max(omega(iMode)-2,0.5);
wEnd = omega(iMode)+2;

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

As = [1 3 10 30];

for i = 1:length(As)
    file = [fileparts(which(mfilename)) filesep 'FRF_' sprintf('A = %0.1f',As(i)) '.mat'];
    if ~isfile(file)
        results = hbm_frf(hbm,problem,As(i),w0,[],wEnd,[]);
        save(file,'-struct','results');
    else
        results = load(file);
    end

    x = results.X;
    u = results.U;
    w = results.w;
    
    X = permute((x(2,:,:)),[2 3 1]);
    U = permute((u(2,:,:)),[2 3 1]);
    H = X(P.iDof,:)./U(P.iInput,:);
    plot(ax(1),w,abs(H));
    plot(ax(2),w,unwrap(angle(H),[],2));
    
    [~,ii] = max(H);
    sol = hbm_res(hbm,problem,w(ii),As(i),x(:,:,ii));
    Xres(:,:,i) = sol.X;
    wres(i) = sol.w0;
    plot(ax(1),sol.w0,abs(sol.H),'o')
    plot(ax(2),sol.w0,angle(sol.H),'o')
    drawnow
end

bb = hbm;
bb.cont.method = 'none';

A0 = As(1);     w0 = wres(1);     X0 = Xres(:,:,1);
Aend = As(end); wEnd = wres(end); XEnd = Xres(:,:,end);
tic;
sol = hbm_bb(bb,problem,A0,w0,X0,Aend,wEnd,XEnd);
toc;
plot(ax(1),sol.w,abs(sol.H));
plot(ax(2),sol.w,unwrap(angle(sol.H),[],2));