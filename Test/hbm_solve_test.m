function hbm_solve_test
problem = test_params;

hbm.harm.Nfft = 32;

hbm.options.bUseStandardHBM = true;
hbm.dependence.x = true;
hbm.dependence.xdot  = true;
hbm.dependence.xddot  = true;
hbm.dependence.u = true;
hbm.dependence.udot  = true;
hbm.dependence.uddot  = true;
hbm.dependence.w = true;

hbm.options.bAnalyticalDerivs = true;
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

omega = sqrt(eig(problem.K,problem.M));
w0 = 4;
A = 100;

% hbm.harm.NHarm = 4;
hbm.harm.group{1}.kHarm = [0; 1];
hbm.harm.group{2}.kHarm = [0; 1; 2; 3];
% 
problem.iGroup = [1 2]';

[hbm,problem] = setuphbm(hbm,problem);
NDof = problem.NDof;

hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';
hbm = setuphbm(hbm,problem);
tic;
sol1 = hbm_solve(hbm,problem,w0,A);
[t1,x1,xdot1] = get_time_series(hbm,w0,sol1.X);
tRun(1) = toc;

hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'sum';
hbm = setuphbm(hbm,problem);
tic;
sol2 = hbm_solve(hbm,problem,w0,A);
[t2,x2,xdot2] = get_time_series(hbm,w0,sol2.X);
tRun(2) = toc;

hbm.options.aft_method = 'fft';
hbm.options.jacob_method = 'mat';
hbm = setuphbm(hbm,problem);
tic;
sol3 = hbm_solve(hbm,problem,w0,A);
[t3,x3,xdot3] = get_time_series(hbm,w0,sol3.X);
tRun(3) = toc;

hbm.options.aft_method = 'fft';
hbm.options.jacob_method = 'sum';
hbm = setuphbm(hbm,problem);
tic;
sol4 = hbm_solve(hbm,problem,w0,A);
[t4,x4,xdot4] = get_time_series(hbm,w0,sol4.X);
tRun(4) = toc;

tic;
y0 = [x4(1,:) xdot4(1,:)];
fun = @(t,y)test_odefun(t,y,w0,sol1.U,hbm,problem);
M = blkdiag(eye(NDof),problem.M);
options = odeset('Mass',M,'Vectorized','on','MassSingular','yes','RelTol',1E-12,'AbsTol',1E-12);
[t5,y] = ode15s(fun,[t4(1) t4(end)],y0,options);
tRun(5) = toc;

x5 = y(:,1:NDof);
xdot5 = y(:,NDof+(1:NDof));

figure
subplot(2,2,1)
plot(t1,x1(:,1:NDof)','-',t2,x2(:,1:NDof)','-',t3,x3(:,1:NDof)','-',t4,x4(:,1:NDof)','-',t5,x5(:,1:NDof)','o-');
ylabel('x');

subplot(2,2,3)
plot(t1,xdot1(:,1:NDof)','-',t2,xdot2(:,1:NDof)','-',t3,xdot3(:,1:NDof)','-',t4,xdot4(:,1:NDof)','-',t5,xdot5(:,1:NDof)','o-');
ylabel('xdot');
xlabel('Time (s)');

names = {'HBM-mat/mat','HBM-mat/sum','HBM-fft/mat','HBM-fft/sum','ODE'};

subplot(1,2,2)
plot(x1(:,1:NDof),xdot1(:,1:NDof));
hold on
plot(x2(:,1:NDof),xdot2(:,1:NDof));
plot(x3(:,1:NDof),xdot3(:,1:NDof));
plot(x4(:,1:NDof),xdot4(:,1:NDof));
plot(x5(:,1:NDof),xdot5(:,1:NDof),'o-');
xlabel('x');
ylabel('xdot');
leg = {};
for i = 1:5
    leg = [leg cellsprintf([names{i} '(x%d)'],num2cell(1:NDof))];
end
legend(leg)
    
for i = 1:5
     fprintf('%s: %f s\n',names{i},tRun(i));
end

sol1.X
sol2.X
sol3.X
sol4.X