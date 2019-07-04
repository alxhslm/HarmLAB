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
if hbm.options.bUseStandardHBM && hbm.harm.NHarm(2) > 0
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
hbm.harm = setupGroups(problem,hbm.harm);

if ~isfield(problem,'sparsity')
    hbm.sparsity = ones(hbm.harm.NComp*problem.NDof);
else
    hbm.sparsity = repmat(problem.sparsity(1:problem.NDof,1:problem.NDof),hbm.harm.NComp);
end
iRetain = hbm.harm.iRetainNL;
hbm.sparsity = hbm.sparsity(iRetain,iRetain);

%% Precompute matrices
hbm.lin    = setupLin(hbm.harm,problem);
hbm.nonlin = setupNonlin(hbm.harm,problem);

hbm.harm = default_missing(hbm.harm,{'iHarmPlot'},{1:hbm.harm.NFreq});

if ~isfield(problem,'iDofPlot') && ~isfield(problem,'RDofPlot')
    problem.iDofPlot = 1:problem.NDof;
end

if ~isfield(problem,'RDofPlot')
    R = zeros(length(problem.iDofPlot),problem.NDof);
    for j = 1:length(problem.iDofPlot)
        R(j,j) = 1;
    end
    problem.RDofPlot = R;
elseif size(problem.RDofPlot,2) ~= problem.NDof
    error('Wrong size for RDofPlot')
end
  

function options = setupOptions(options)
if ~isfield(options,'bAnalyticalDerivs'), options.bAnalyticalDerivs = 1; end
if ~isfield(options,'bUseStandardHBM'), options.bUseStandardHBM = 0; end
if ~isfield(options,'solver'), options.solver = 'ipopt'; end
options = default_missing(options,{'aft_method','jacob_method'},{'mat','mat'});

function dependence = setupDependence(dependence)
dependence = default_missing(dependence,{'x','xdot','xddot','w','u','udot','uddot'},{true,false,false,false,false,false,false});

function scaling = setupScaling(scaling)
scaling = default_missing(scaling,{'method','tol'},{'max',1E-6});


function cont = setupCont(cont)
cont = default_missing(cont,{'method','bUpdate','step0','min_step','max_step','ftol','xtol','c', 'C','maxfail','num_iter_increase','num_iter_reduce'},{'predcorr',true,1E-3, 1E-6, 5E-3, 1E-6,1E-6,0.5, 1.05,4,10,3});

if ~isfield(cont,'predcorr'), cont.predcorr = struct(); end
cont.predcorr = default_missing(cont.predcorr,{'predictor','corrector','bMoorePenrose','solver','maxit'},{'linear','pseudo',true,'ipopt',30});

if ~isfield(cont,'coco'), cont.coco = struct(); end
cont.coco = default_missing(cont.coco,{'ItMX','NPR'},{2E4,5});


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
    f = {'input','output','iInput','iOutput','iHarm'};

    if ~isfield(problem.res,f{i})
        error('Missing field %s from resonance condition',f{i})
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

function harm = setupGroups(problem,harm)
%asign each DOF to a group
for i = 1:length(harm.group)
    harm.group{i}.iDof = find(problem.iGroup==i);
    harm.group{i}.NDof = length(harm.group{i}.iDof);
    harm.group{i}.NCompTot = harm.group{i}.NDof * harm.group{i}.NComp;
end

%now we need to work out which indices we need to retain
NDof = problem.NDof;
NNL  = problem.NNL;
bRetain = false(harm.NComp*NDof,1);
bRetainNL = false(harm.NComp*NNL,1);
for i = 1:length(harm.group)
    for k = 1:harm.group{i}.NFreq
        iDof = harm.group{i}.iDof;
        iNL = dof2nl(problem.iNL,iDof);
        if harm.group{i}.iFreq(k) == 1
            iKeep = iDof;
            iKeepNL = iNL;
        else
            iKeep = NDof + (harm.group{i}.iFreq(k)-2)*2*NDof + [iDof; NDof + iDof];
            iKeepNL = NNL + (harm.group{i}.iFreq(k)-2)*2*NNL + [iNL; NNL + iNL];
        end
        bRetain(iKeep) = true;
        bRetainNL(iKeepNL) = true;
    end
end
%index into the FULL phasor vector comprising of:
%  - including ALL DOFs (linear and non-linear)
%  - at ALL possible harmonics (even ones being neglected for certain DOF)
harm.iRetain = find(bRetain);
harm.NRetain = sum(bRetain);

%index into the PARTIAL phasor vector comprising of:
%  - only the non-linear DOFs
%  - at ALL possible harmonics (even ones being neglected for certain DOF)
harm.iRetainNL = find(bRetainNL);
harm.NRetainNL = sum(bRetainNL);

bLin = false(problem.NDof*harm.NComp,1);
bNL = false(problem.NDof*harm.NComp,1);
for j = 1:harm.NFreq
    if j == 1
        iLin  = problem.iLin;
        iNL  = problem.iNL;
    else
        iLin = NDof + 2*(j-2)*problem.NDof + [problem.iLin; NDof + problem.iLin];
        iNL  = NDof + 2*(j-2)*problem.NDof + [problem.iNL; NDof + problem.iNL];
    end
    bLin(iLin) = true;
    bNL(iNL) = true;
end
harm.iLin = find(bLin(bRetain));
harm.iNL = find(bNL(bRetain));

function j = dof2nl(nonlin,ind)
j = [];
for i = 1:length(ind)
    j = [j; find(nonlin == ind(i))];
end

function s = default_missing(s,f,d)
for i = 1:length(f)
    if ~isfield(s,f{i})
        s.(f{i}) = d{i};
    end
end