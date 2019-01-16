function problem = hbm_scaling(problem,hbm,x,w,A)
X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
tol = hbm.scaling.tol;
switch hbm.scaling.method
    case 'max'
        xdc  = max(abs(X(1,:)),tol);
        xmax = max(max(abs(X(2:end,:)),[],1),tol);
        xharm = repmat(xmax,hbm.harm.NFreq-1,1);
    case 'abs'
        xdc = max(abs(X(1,:)),tol);
        xharm = max(abs(X(2:end,:)),tol);
end
xscale = [xdc; xharm*(1+1i)];
problem.xscale = packdof(xscale,hbm.harm.iRetain);
problem.Xscale = problem.xscale;
problem.Fscale = ones(length(problem.xscale),1);

if nargin > 3
    %frf,resonance, or bb
    problem.wscale = w;
    problem.Xscale(end+1) = w;
end

if nargin > 4
    %bb
    problem.Ascale = A;
    problem.Xscale(end+1) = A;
    
    %need extra resonance condition
    problem.Fscale(end+1) = 1;
end

problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';

% %from solve
% X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
% tol = 1E-6;
% xdc  = max(abs(X(1,:)),tol);
% xmax = max(max(abs(X(2:end,:)),[],1),tol);
% xharm = repmat(xmax,hbm.harm.NFreq-1,1);
% % xharm = max(abs(X(2:end,:)),1E-10);
% xscale = [xdc; xharm*(1+1i)];
% problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
% problem.Fscale = x*0+1;
% problem.Jscale = (1./problem.Fscale(:))*problem.xscale(:)';
% 
% %from frf
% X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
% xdc  = max(abs(X(1,:)),1E-6);
% xmax = max(max(abs(X(2:end,:)),[],1),1E-6);
% xharm = repmat(xmax,hbm.harm.NFreq-1,1);
% % xharm = max(abs(X(2:end,:)),1E-10);
% xscale = [xdc; xharm*(1+1i)];
% problem.xscale = packdof(xscale,hbm.harm.iRetain);
% problem.wscale = w;
% problem.Fscale = x*0+1;
% problem.Xscale = [problem.xscale; problem.wscale];
% problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';
% 
% %from resonnace
% X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
% xdc  = max(abs(X(1,:)),1E-6);
% xharm = repmat(max(abs(X(2:end,:)),[],1),hbm.harm.NFreq-1,1);
% xharm = max(xharm,1E-6);
% xscale = [xdc; xharm*(1+1i)];
% problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
% problem.Fscale = x*0+1;
% problem.wscale = w0;
% problem.Xscale = [problem.xscale; problem.wscale];
% problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';
% 
% %from bb
% X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
% xdc  = max(abs(X(1,:)),1E-6);
% xmax = max(max(abs(X(2:end,:)),[],1),1E-6);
% xharm = repmat(xmax,hbm.harm.NFreq-1,1);
% % xharm = max(abs(X(2:end,:)),1E-10);
% xscale = [xdc; xharm*(1+1i)];
% problem.xscale = packdof(xscale,hbm.harm.iRetain);
% problem.wscale = 1;%w;
% problem.Ascale = 1E-4;%A;
% problem.Fscale = ones(length(problem.xscale)+problem.constr.N,1);
% problem.Xscale = [problem.xscale; problem.wscale; problem.Ascale];
% problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';
