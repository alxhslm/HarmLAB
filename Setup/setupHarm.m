function harm = setupHarm(harm)
%sort out order of harmonics

if ~isfield(harm,'group')
    f = {'NHarm','kHarm'};
    for i = 1:length(f)
        if isfield(harm,f{i})
            harm.group{1}.(f{i}) = harm.(f{i});
        end
    end
end 

if ~iscell(harm.group)
    harm.group = {harm.group};
end

for i = 1:length(harm.group)
    [harm.group{i},NHarmGroup(i,:)] = setupGroup(harm.group{i});
end
harm.NHarm = max(NHarmGroup,[],1);

%now default missing fields
f = {'Nfft','rFreqRatio','rFreqBase'};
d = {max(8*harm.NHarm,1),0*harm.NHarm+1,0*harm.NHarm+1};
for i = 1:length(f)
    if ~isfield(harm,f{i})
        harm.(f{i}) = d{i};
        warning('Field %s is missing from the options. Defaulting to %s',f{i},mat2str(d{i}))
    end
end

%check all the sizes match
if length(harm.rFreqRatio) ~= length(harm.NHarm)
    error('The number of frequencies should match the number of elements in the NHarm vector')
end

if length(harm.rFreqRatio) ~= length(harm.Nfft)
    error('The number of frequencies should match the number of elements in the Nfft vector')
end

%now add dummy second harmonic for single harmonic problems
if length(harm.rFreqRatio) == 1
    harm.rFreqRatio(2) = 1;
    harm.NHarm(2) = 0;
    harm.Nfft(2) = 1;
    for i = 1:length(harm.group)
        harm.group{i}.kHarm(:,2) = 0;
    end
end

%deal with any sub-harmonics, changing base frequencies accordingly
harm = detect_subharmonics(harm);

harm.kHarm = zeros(0,2);
for i = 1:length(harm.group)
    %put 0th harmonic first
    harm.group{i}.kHarm = reorder_harmonics(harm.group{i}.kHarm);
    
    harm.group{i}.NFreq = size(harm.group{i}.kHarm,1);
    harm.group{i}.NComp = (2*(harm.group{i}.NFreq-1)+1);
    
    harm.kHarm = [harm.kHarm;
                  harm.group{i}.kHarm];
end

harm.kHarm = reorder_harmonics(unique(harm.kHarm,'rows','stable'));
 %find the subscripts for using FFT/IFFT
harm.iSub  = find_subscripts(harm.kHarm,harm.NHarm);
harm.NFreq = size(harm.kHarm,1);
harm.NComp = (2*(harm.NFreq-1)+1);

for i = 1:length(harm.group)
    for k = 1:harm.group{i}.NFreq
        harm.group{i}.iFreq(k) = find(harm.kHarm(:,1)==harm.group{i}.kHarm(k,1) & harm.kHarm(:,2)==harm.group{i}.kHarm(k,2));
    end
end

function [group,NHarm] = setupGroup(group)
if isfield(group,'kHarm')
    group.NHarm = max(abs(group.kHarm),[],1);
elseif isfield(group,'NHarm')
    if length(group.NHarm) == 1
        group.kHarm = (0:group.NHarm)';
    else
        k1 = -group.NHarm(1):group.NHarm(1);
        k2 = -group.NHarm(2):group.NHarm(2);
        [K1,K2] = ndgrid(k1,k2);
        keep = (K1+K2)>0 | ((K1+K2)== 0 & K1>=0);
        group.kHarm = [K1(keep)  K2(keep)];
    end
else
    error('Need either kHarm or NHarm');
end
NHarm = max(abs(group.kHarm),[],1);
group = rmfield(group,'NHarm');

function harm = detect_subharmonics(harm)
ratio = harm.NHarm*0 + 1;

for i = 1:length(harm.group)
    [~,fac] = rat(harm.group{i}.kHarm);

    for k = 1:2
        for j = 1:size(harm.group{i}.kHarm,1)
            ratio(k) = lcm(fac(j,k),ratio(k));
        end
    end
end
for i = 1:length(harm.group)
    harm.group{i}.kHarm = harm.group{i}.kHarm .* repmat(ratio,size(harm.group{i}.kHarm,1),1);
end
harm.rFreqBase = harm.rFreqBase./ratio;
harm.NHarm = harm.NHarm.*ratio;
harm.Nfft = harm.Nfft.*ratio;

function kHarm = reorder_harmonics(kHarm)
%find 0th harmonic and put it first

zero  = kHarm(:,1) == 0 & kHarm(:,2) == 0;
first1 = kHarm(:,1) == 1 & kHarm(:,2) == 0;
first2 = kHarm(:,1) == 0 & kHarm(:,2) == 1;
rest = ~zero & ~first1 & ~first2;
ii = [find(zero); find(first1); find(first2); find(rest)];
kHarm = kHarm(ii,:);

function iSub = find_subscripts(kHarm,NHarm)
k1 = -NHarm(1):NHarm(1);
k2 = -NHarm(2):NHarm(2);
[K1,K2] = ndgrid(k1,k2);
iSub = 0*kHarm;
for i = 1:size(kHarm,1)
    [iSub(i,1),iSub(i,2)] = ind2sub([2*NHarm(1)+1,2*NHarm(2)+1],find(K1 == kHarm(i,1) & K2 == kHarm(i,2)));
end