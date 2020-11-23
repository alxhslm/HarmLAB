function hbm_bb_test(b3d)
problem = test_params;

P = problem.P;

hbm.harm.rFreqRatio = 1;
hbm.harm.NHarm = 2;
hbm.harm.Nfft = 32;
hbm.harm.iHarmPlot = 2;

hbm.dependence.x = true;
hbm.dependence.xdot  = true;
hbm.dependence.w = false;

hbm.scaling.tol = 1;

hbm.cont.step0 = 1E-2;
hbm.cont.max_step = 1E-1;
hbm.cont.min_step = 1E-6;

if nargin < 1
    b3d = 1;
end

if b3d
    hbm.harm.rFreqRatio(end+1) = 1.376;
    hbm.harm.NHarm(end+1) = 2;
    hbm.harm.Nfft(end+1) = 16;
    hbm.harm.iHarmPlot(end+1) = 3;
end

hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'arclength';

problem.res.iOutput = P.iDof;
problem.res.iInput = P.iInput;
problem.res.iHarm = 2;
problem.res.sign = 1;
problem.res.output = 'x';
problem.res.input = 'u';

[hbm,problem] = setuphbm(hbm,problem);

iRes = [];
iRes(end+1) = find(hbm.harm.kHarm(:,1) == 1 & hbm.harm.kHarm(:,2) == 0);
iRes(end+1) = find(hbm.harm.kHarm(:,1) == 0 & hbm.harm.kHarm(:,2) == 1);
NRes = length(iRes);

iMode = 1;
omega = sqrt(eig(problem.K,problem.M));
w0 = max(omega(iMode)-2,0.5);
wEnd = omega(iMode)+2;

figure
for j = 1:NRes
    ax(1,j) = subplot(2,NRes,j);
    hold on
    if j == 1
        ylabel('|X| (mag)');
    end
    title(sprintf('iHarm %d',iRes(j)))
    
    ax(2,j) = subplot(2,NRes,NRes+j);
    hold on
    if j == 1
        ylabel('\angle X (deg)');
    end
    
    xlabel('\omega (rads)');
end
linkaxes(ax(:),'x')
xlim([w0,wEnd]);

As = [1 2 3 4];

root = fileparts(which(mfilename));
for i = 1:length(As)
    file = fullfile(root,sprintf('FRF_A = %0.1f.mat',As(i)));
    if ~isfile(file)
        sol = hbm_frf(hbm,problem,As(i),w0,[],wEnd,[]);
        results.X = cat(3,sol.X);
        results.U = cat(3,sol.U);
        results.w = cat(2,sol.w);
        save(file,'-struct','results');
    else
        results = load(file);
    end

    x = results.X;
    u = results.U;
    w = results.w;
    
    for j = 1:NRes
        X = permute((x(iRes(j),:,:)),[2 3 1]);
        U = permute((u(iRes(j),:,:)),[2 3 1]);
        H = X(P.iDof,:)./U(P.iInput,:);
        plot(ax(1,j),w,abs(H));
        plot(ax(2,j),w,unwrap(angle(H),[],2));
    
        file = fullfile(root,sprintf('RES_A = %0.1f_%d.mat',As(i),j));
        if ~isfile(file)
            problem.res.iHarm = iRes(j);
            [hbm,problem] = setuphbm(hbm,problem);
            [~,ii] = max(H);
            res = hbm_res(hbm,problem,w(ii),As(i),x(:,:,ii));
            save(file,'-struct','res');
        else
            res = load(file);
        end

        Xres{j}(:,:,i) = res.X;
        wres{j}(i) = res.w;
        plot(ax(1,j),res.w,abs(res.H),'o')
        plot(ax(2,j),res.w,angle(res.H),'o')
    end
    drawnow
end

bb = hbm;
bb.cont.method = 'none';

for j = 1:NRes
    A0 = As(1);     w0 = wres{j}(1);     X0 = Xres{j}(:,:,1);
    Aend = As(end); wEnd = wres{j}(end); XEnd = Xres{j}(:,:,end);
    
    problem.res.iHarm = iRes(j);
    [hbm,problem] = setuphbm(hbm,problem);
    
    file = fullfile(root,sprintf('BB_%d.mat',j));
    if ~isfile(file)
        sol = hbm_bb(bb,problem,A0,w0,X0,Aend,wEnd,XEnd);
        results.X = cat(3,sol.X);
        results.U = cat(3,sol.U);
        results.w = cat(2,sol.w);
        results.H = cat(2,sol.H);
        save(file,'-struct','results');
    else
        results = load(file);
    end
    
    for k = 1:NRes
        X = squeeze(results.X(iRes(k),P.iDof,:));
        U = squeeze(results.U(iRes(k),P.iInput,:));
        H = X./U;
        plot(ax(1,k),results.w,abs(H));
        plot(ax(2,k),results.w,unwrap(angle(H)));
    end
end