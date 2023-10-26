function hbm_solve_test(b3d)
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

omega = sqrt(eig(problem.K,problem.M));
w0 = 4;
A = 100;

[hbm,problem] = setuphbm(hbm,problem);
NDof = problem.NDof;

tic;
sol1 = hbm_solve(hbm,problem,w0,A,[]);
[t1,x1,xdot1] = get_time_series3d(hbm,w0,sol1.X);
tRun(1) = toc;

tic;
y0 = [x1(1,:) xdot1(1,:)];
fun = @(t,y)test_odefun(t,y,w0*hbm.harm.rFreqRatio,sol1.U,hbm,problem);
M = blkdiag(eye(NDof),problem.M);
options = odeset('Mass',M,'Vectorized','on','MassSingular','yes','RelTol',1E-12,'AbsTol',1E-12);
[t2,y] = ode15s(fun,[t1(1) t1(end)],y0,options);
tRun(2) = toc;

x2 = y(:,1:NDof);
xdot2 = y(:,NDof+(1:NDof));

figure
subplot(2,2,1)
plot(t1,x1(:,1:NDof)','-',t2,x2(:,1:NDof)','o-');
ylabel('x');

subplot(2,2,3)
plot(t1,xdot1(:,1:NDof)','-',t2,xdot2(:,1:NDof)','o-');
ylabel('xdot');
xlabel('Time (s)');

names = {'HBM','ODE'};

subplot(1,2,2)
plot(x1(:,1:NDof),xdot1(:,1:NDof));
hold on
plot(x2(:,1:NDof),xdot2(:,1:NDof),'o-');
xlabel('x');
ylabel('xdot');
leg = {};
for i = 1:length(tRun)
    leg = [leg cellsprintf([names{i} '(x%d)'],num2cell(1:NDof))];
end
legend(leg)
    
for i = 1:length(tRun)
     fprintf('%s: %f s\n',names{i},tRun(i));
end