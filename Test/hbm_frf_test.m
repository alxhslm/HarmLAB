function hbm_frf_test(b3d)
problem = test_params;

hbm.harm.rFreqRatio = 1;
hbm.harm.NHarm = 2;
hbm.harm.Nfft = 32;
hbm.harm.iHarmPlot = 2;

if nargin < 1
    b3d = 1;
end

if b3d
    hbm.harm.rFreqRatio(end+1) = 1.376;
    hbm.harm.NHarm(end+1) = 2;
    hbm.harm.Nfft(end+1) = 16;
    hbm.harm.iHarmPlot(end+1) = 3;
end

hbm.dependence.x = true;
hbm.dependence.xdot  = true;
hbm.dependence.w = false;

hbm.scaling.tol = 1;

hbm.cont.step0 = 1E-2;
hbm.cont.max_step = 1E-1;
hbm.cont.min_step = 1E-6;

hbm.cont.method = 'predcorr';

[hbm,problem] = setuphbm(hbm,problem);

omega = sqrt(eig(problem.K,problem.M));
w0 = max(min(omega)-2,0.5);
wEnd = max(omega)+2;
A = 1;

S = {};

tic;
hbm.cont.method = 'none';
sol = hbm_frf(hbm,problem,A,w0,[],wEnd,[]);
S{end+1} = storeResults(sol,toc,'none');

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'arclength';
sol = hbm_frf(hbm,problem,A,w0,[],wEnd,[]);
S{end+1} = storeResults(sol,toc,'arclength');

tic;
hbm.cont.method = 'predcorr';
hbm.cont.predcorr.corrector = 'pseudo';
sol = hbm_frf(hbm,problem,A,w0,[],wEnd,[]);
S{end+1} = storeResults(sol,toc,'pseudo');

tic;
S{end+1} = ode_frf(w0,wEnd,A,hbm,problem);

figure
hold on
for j = 1:problem.NDof
    ax_mag(j) = subplot(2,problem.NDof,j);
    hold on
    for i = 1:length(S)
        hmag(i,j) = plot(ax_mag(j),S{i}.w,abs(S{i}.X(j,:)));
    end
    ylabel(ax_mag(j),sprintf('|X_%d| (mag)',j));
end

for j = 1:problem.NDof
    ax_ph(j) = subplot(2,problem.NDof,problem.NDof+j);
    hold on
    for i = 1:length(S)
        hph(i,j) = plot(ax_ph(j),S{i}.w,unwrap(angle(S{i}.X(j,:)),[],2));
    end
    xlabel(ax_ph(j),'\omega (rads)');
    ylabel(ax_ph(j),sprintf('\\angle X_%d (deg)',j));
end
linkaxes([ax_mag ax_ph],'x')

for i = 1:length(S)
    leg{i} = S{i}.name;
end
legend(ax_ph(end),hph(:,end),leg)

for i = 1:length(S)
    fprintf('%10s : %0.2f s\n',S{i}.name,S{i}.t)
end

function S = ode_frf(w0,wEnd,A,hbm,problem)
y = zeros(1,2*problem.NDof);
ws = linspace(w0,wEnd,20);
NCycle = 50;
Nfft = 101;
for i = 1:length(ws)
    y0 = y(end,:);
    w0 = hbm.harm.rFreqRatio * ws(i);
    U = A*test_excite(hbm,problem,w0);
    T = 2*pi/ws(i);
    [t,y] = ode45(@(t,y)test_odefun(t,y,w0,U,hbm,problem),[0 NCycle*T],y0);
    
    ii = find(t>(t(end)-T),1);
    t = t(ii:end) - t(ii);
    y = y(ii:end,:);
    
    tReg = linspace(0,T,Nfft)';
    Fs = Nfft/T;
    yReg = interp1(t,y,tReg,'linear','extrap');
    
    Yfft = 2*fft(yReg,[],1)/Nfft;
    Wfft = (0:(Nfft-1))'/Nfft * 2*pi* Fs;
        
    X(:,i) = Yfft(2,1:problem.NDof);
    w(i) = Wfft(2);
end

S.X = X;
S.w = w;
S.t = toc;
S.name = 'ode';

function S = storeResults(sol,t,name)
X = cat(3,sol.X);

S.X = permute(X(2,:,:),[2 3 1]);
S.w = cat(2,sol.w);
S.t = t;
S.name = name;