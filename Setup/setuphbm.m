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
problem = setupProblem(problem);
hbm.harm = setupGroups(problem,hbm.harm);

if ~isfield(problem,'sparsity')
    hbm.sparsity = ones(hbm.harm.NComp*problem.NDof + prod(hbm.harm.Nfft)*problem.NAlg);
else
    Jxx = repmat(problem.sparsity(1:problem.NDof,1:problem.NDof),hbm.harm.NComp);
    Jxa = repmat(problem.sparsity(1:problem.NDof,problem.NDof+1:end),hbm.harm.NComp,prod(hbm.harm.Nfft));
    Jax = repmat(problem.sparsity(problem.NDof+1:end,1:problem.NDof),prod(hbm.harm.Nfft),hbm.harm.NComp);
    Jaa = blkmat(repmat(problem.sparsity(problem.NDof+1:end,problem.NDof+1:end),1,1,prod(hbm.harm.Nfft)));

    iRetain = hbm.harm.iRetain;
    hbm.sparsity = [Jxx(iRetain,iRetain) Jxa(iRetain,:);
                    Jax(:,iRetain) Jaa];
end

if problem.NAlg > 0 && strcmp(hbm.options.jacob_method,'sum')
    warning('Defaulting Jacobian method to "mat". "sum" is not supported for problems with algebraic constraints')
    hbm.options.jacob_method = 'mat';
end

%% Precompute matrices
hbm.lin    = setupLin(hbm.harm,problem);
hbm.nonlin = setupNonlin(hbm.harm,problem);

hbm.harm = default_missing(hbm.harm,{'iHarmPlot'},{1:hbm.harm.NFreq});
problem = default_missing(problem,{'iDofPlot'},{1:problem.NDof});

function options = setupOptions(options)
if ~isfield(options,'bAnalyticalDerivs'), options.bAnalyticalDerivs = 1; end
if ~isfield(options,'bUseStandardHBM'), options.bUseStandardHBM = 0; end
if ~isfield(options,'solver'), options.solver = 'ipopt'; end
options = default_missing(options,{'aft_method','jacob_method'},{'mat','mat'});

function dependence = setupDependence(dependence)
dependence = default_missing(dependence,{'x','xdot','xddot','w','u','udot','uddot','xalg'},{true,false,false,false,false,false,false,false});

function scaling = setupScaling(scaling)
scaling = default_missing(scaling,{'method','tol'},{'max',1E-6});


function cont = setupCont(cont)
cont = default_missing(cont,{'method','bUpdate','step0','min_step','max_step','ftol','xtol','c', 'C','maxfail','num_iter_increase','num_iter_reduce'},{'predcorr',true,1E-3, 1E-6, 5E-3, 1E-6,1E-6,0.5, 1.05,4,10,3});

if ~isfield(cont,'predcorr'), cont.predcorr = struct(); end
cont.predcorr = default_missing(cont.predcorr,{'predictor','corrector','bMoorePenrose','solver','maxit'},{'linear','pseudo',true,'ipopt',30});

if ~isfield(cont,'coco'), cont.coco = struct(); end
cont.coco = default_missing(cont.coco,{'ItMX','NPR'},{2E4,5});


function problem = setupProblem(problem)
if ~isfield(problem,'name')
    problem.name = '';
end
if ~isfield(problem,'NDof')
    problem.NDof = size(problem.K,2);
end
if ~isfield(problem,'NAlg')
    problem.NAlg = 0;
end
if ~isfield(problem,'NInput')
    problem.NInput = size(problem.Ku,2);
end
if ~isfield(problem,'NOutput')
    error('NOutput is missing from problem structure')
end

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

if ~isfield(problem,'iGroup')
    problem.iGroup = ones(problem.NDof,1);
end

function harm = setupGroups(problem,harm)
%asign each DOF to a group
for i = 1:length(harm.group)
    harm.group{i}.iDof = find(problem.iGroup==i);
    harm.group{i}.NDof = length(harm.group{i}.iDof);
    harm.group{i}.NCompTot = harm.group{i}.NDof * harm.group{i}.NComp;
end

%now we need to work out which indices we need to retain
NDof = problem.NDof;
bRetain = zeros(harm.NComp*NDof,1);
for i = 1:length(harm.group)
    for k = 1:harm.group{i}.NFreq
        iDof = harm.group{i}.iDof;
        if harm.group{i}.iFreq(k) == 1
            iKeep = iDof;
        else
            iKeep = NDof + (harm.group{i}.iFreq(k)-2)*2*NDof + [iDof; NDof + iDof];
        end
        bRetain(iKeep) = 1;
    end
end
harm.iRetain = find(bRetain);
harm.NRetain = sum(bRetain);

function s = default_missing(s,f,d)
for i = 1:length(f)
    if ~isfield(s,f{i})
        s.(f{i}) = d{i};
    end
end