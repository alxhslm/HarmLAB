function [hbm,problem] = setuphbm(hbm,problem)

%% Harmonics
if ~isfield(hbm,'harm')
    hbm.harm = struct();
end
hbm.harm = setupHarm(hbm.harm);

%% Options
if ~isfield(hbm,'options')
    hbm.options = struct();
end
hbm.options = setupOptions(hbm.options);
if hbm.options.bUseStandardHBM && prod(hbm.harm.NHarm) > 0
    error('Cannot use standard HBM code in case of more than one fundemental')
end

%% Dependence
if ~isfield(hbm,'dependence')
    hbm.dependence = struct();
end
hbm.dependence = setupDependence(hbm.dependence);

%% Scaling
if ~isfield(hbm,'scaling')
    hbm.scaling = struct();
end
hbm.scaling = setupScaling(hbm.scaling);

%% Continuation
if ~isfield(hbm,'cont')
    hbm.cont = struct(); 
end
hbm.cont = setupCont(hbm.cont);

%% Problem definition
problem = setupProblem(problem,hbm);
hbm.harm = setupNL(problem,hbm.harm);

if ~isfield(problem,'sparsity')
    hbm.sparsity = ones(hbm.harm.NComp*problem.NNL);
else
    hbm.sparsity = repmat(problem.sparsity(problem.iNL,problem.iNL),hbm.harm.NComp);
end

%% Precompute matrices
hbm.lin    = setupLin(hbm.harm,problem);
hbm.nonlin = setupNonlin(hbm.harm,problem);

hbm.harm = default_missing(hbm.harm,{'iHarmPlot'},{1:hbm.harm.NFreq});

if ~isfield(problem,'RDofPlot')
    if ~isfield(problem,'iDofPlot')
        problem.iDofPlot = 1:problem.NDof;
    end

    R = zeros(length(problem.iDofPlot),problem.NDof);
    for j = 1:length(problem.iDofPlot)
        R(j,problem.iDofPlot(j)) = 1;
    end
    
    problem.RDofPlot = R;
elseif size(problem.RDofPlot,2) ~= problem.NDof
    error('Wrong size for RDofPlot')
end
  
function options = setupOptions(options)
if ~isfield(options,'bAnalyticalDerivs'), options.bAnalyticalDerivs = 1; end
if ~isfield(options,'bUseStandardHBM'), options.bUseStandardHBM = 0; end
if ~isfield(options,'solver'), options.solver = 'ipopt'; end

function dependence = setupDependence(dependence)
dependence = default_missing(dependence,{'x','xdot','xddot','w','u','udot','uddot'},{true,false,false,false,false,false,false});

function scaling = setupScaling(scaling)
scaling = default_missing(scaling,{'method','tol'},{'max',1E-6});

function cont = setupCont(cont)
cont = default_missing(cont,{'method','bUpdate','step0','min_step','max_step','ftol','xtol','c', 'C','maxfail','num_iter_increase','num_iter_reduce'},{'predcorr',true,1E-3, 1E-6, 5E-3, 1E-6,1E-6,0.5, 1.05,4,10,3});

if ~isfield(cont,'predcorr'), cont.predcorr = struct(); end
cont.predcorr = default_missing(cont.predcorr,{'predictor','corrector','bMoorePenrose','solver','maxit'},{'linear','pseudo',true,'ipopt',30});

function problem = setupProblem(problem,hbm)
if ~isfield(problem,'name')
    problem.name = '';
end
problem.NDof = size(problem.K,2);
problem.NInput = size(problem.Ku,2);

if ~isfield(problem,'iNL')
    problem.iNL = (1:problem.NDof)';
end
problem.NNL = length(problem.iNL);

tmp = true(problem.NDof,1);
tmp(problem.iNL) = false;
problem.iLin = find(tmp);
problem.NLin = length(problem.iLin);

f = {'K','M','C'};
for i = 1:length(f)
    if ~isfield(problem,f{i})
        problem.(f{i}) = zeros(problem.NDof);
    elseif size(problem.(f{i}),1) ~= problem.NDof || size(problem.(f{i}),2) ~= problem.NDof 
        error('Wrong size for %s matrix',f{i})
    end
end

if ~isfield(problem,'F0')
    problem.('F0') = zeros(problem.NDof,1);
elseif length(problem.('F0')) ~= problem.NDof
    error('Wrong size for %s matrix','F0')
end


f = {'Ku','Cu','Mu'};
for i = 1:length(f)
    if ~isfield(problem,f{i})
        problem.(f{i}) = zeros(problem.NDof,problem.NInput);
    elseif  size(problem.(f{i}),1) ~= problem.NDof || size(problem.(f{i}),2) ~= problem.NInput 
        error('Wrong size for %s matrix',f{i})
    end
end

try
    States = empty_states(problem);
    out = feval(problem.model,'output',States,hbm,problem);
    problem.NOutput = length(out);
catch
     error('Error detected in non-linear function')
end

if ~isfield(problem,'iGroup')
    problem.iGroup = ones(problem.NDof,1);
end

if isfield(problem,'res')
    f = {'input','output','iHarm'};

    if ~isfield(problem.res,f{i})
        error('Missing field %s from resonance condition',f{i})
    end
    
    switch problem.res.input
        case 'unity'
            problem.res.NInput = 1;
        case 'fe'
            problem.res.NInput = problem.NDof;
        otherwise
           problem.res.NInput = problem.NInput;
    end

    if ~isfield(problem.res,'RInput')
        if ~isfield(problem.res,'iInput')
            problem.res.iInput = 1:problem.res.NInput;
        end
        
        R = zeros(length(problem.res.iInput),problem.res.NInput);
        for j = 1:length(problem.res.iInput)
            R(j,problem.res.iInput(j)) = 1;
        end
        
        problem.res.RInput = R;
    elseif size(problem.res.RInput,2) ~= problem.res.NInput
        error('Wrong size for RInput')
    end
    
    switch problem.res.output
        case 'none'
            problem.res.NOutput = 1;
        otherwise
           problem.res.NOutput = problem.NDof;
    end
    
    if ~isfield(problem.res,'ROutput')
        if ~isfield(problem.res,'iOutput')
            problem.res.iOutput = 1:problem.res.NOutput;
        end
        
        R = zeros(length(problem.res.iOutput),problem.res.NOutput);
        for j = 1:length(problem.res.iOutput)
            R(j,problem.res.iOutput(j)) = 1;
        end
        
        problem.res.ROutput = R;
    elseif size(problem.res.ROutput,2) ~= problem.res.NOutput
        error('Wrong size for ROutput')
    end
    
    
end

function States = empty_states(problem)
States.w0 = NaN;
States.wBase = NaN;
States.t = 0;

States.x = zeros(problem.NDof,1);
States.xdot = zeros(problem.NDof,1);
States.xddot = zeros(problem.NDof,1);

States.u = zeros(problem.NInput,1);
States.udot = zeros(problem.NInput,1);
States.uddot = zeros(problem.NInput,1);

function harm = setupNL(problem,harm)
harm.bNL  = false(problem.NDof*harm.NComp,1);
for j = 1:harm.NFreq
    if j == 1
        iNL  = problem.iNL;
    else
        iNL  = problem.NDof + 2*(j-2)*problem.NDof + [problem.iNL; problem.NDof + problem.iNL];
    end
    harm.bNL(iNL) = true;
end
harm.iNL  = find(harm.bNL);
harm.iLin = find(~harm.bNL);

function s = default_missing(s,f,d)
for i = 1:length(f)
    if ~isfield(s,f{i})
        s.(f{i}) = d{i};
    end
end