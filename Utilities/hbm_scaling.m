function problem = hbm_scaling(problem,hbm,results)
X = results.X;
w = results.w;
A = results.A;

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
Xscale = [xdc; xharm*(1+1i)];
problem.Xscale = packdof(Xscale);
problem.Xnlscale = packdof(Xscale(:,problem.iNL));

problem.Yscale = problem.Xscale;
problem.Zscale = problem.Xnlscale;
problem.Fscale = ones(length(problem.Xnlscale),1);

switch problem.type
    case 'frf'
        problem.wscale = w;
        problem.Yscale(end+1) = w;
        problem.Zscale(end+1) = w;
    case 'resonance'
        problem.wscale = w;
        problem.Yscale(end+1) = w;
        problem.Zscale(end+1) = w;
    case 'bb'
        problem.wscale = w;
        problem.Yscale(end+1) = w;
        problem.Zscale(end+1) = w;
        
        %amplitude is the continuation variable
        problem.Ascale = A;
        problem.Yscale(end+1) = A;
        problem.Zscale(end+1) = A;
        %need extra constraint
        problem.Fscale(end+1) = 1;
    case 'amp'
        problem.Ascale = A;
        problem.Yscale(end+1) = A;
        problem.Zscale(end+1) = A;
    case 'solve'
end
   
problem.Jscale = (1./problem.Fscale(:))*problem.Zscale(:)';