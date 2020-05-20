function harm = setupHarm(harm)

harm = setupHarmonics(harm);

%default missing fields
f = {'Nfft','rFreqRatio','rFreqBase','wFreq0'};
d = {max(8*harm.NHarm,1),0*harm.NHarm+1,0*harm.NHarm+1,0*harm.NHarm};
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
    harm.kHarm(:,2) = 0;
end

%deal with any sub-harmonics, changing base frequencies accordingly
harm = detect_subharmonics(harm);

harm.kHarm = reorder_harmonics(unique(harm.kHarm,'rows','stable'));

%find the subscripts for using FFT/IFFT
harm.NFreq = size(harm.kHarm,1);
harm.NComp = (2*(harm.NFreq-1)+1);

function group = setupHarmonics(group)
if isfield(group,'kHarm')
    group.NHarm = max(abs(group.kHarm),[],1);
elseif isfield(group,'NHarm')
    if length(group.NHarm) == 1
        group.kHarm = (0:group.NHarm)';
    elseif group.NHarm(1) == 0
        k2 = 0:group.NHarm(2);
        k1 = 0*k2;
        group.kHarm = [k1(:) k2(:)];
    elseif group.NHarm(2) == 0
        k1 = 0:group.NHarm(1);
        k2 = 0*k1;
        group.kHarm = [k1(:) k2(:)];
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

function harm = detect_subharmonics(harm)
ratio = harm.NHarm*0 + 1;

[~,fac] = rat(harm.kHarm);

for k = 1:2
    for j = 1:size(harm.kHarm,1)
        ratio(k) = lcm(fac(j,k),ratio(k));
    end
end

harm.kHarm = harm.kHarm .* repmat(ratio,size(harm.kHarm,1),1);

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