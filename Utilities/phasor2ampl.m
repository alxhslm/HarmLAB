function [amp,x] = phasor2ampl(X,Nfft)    
if ~iscell(X)
    X = {X};
end

if nargin < 2
    Nfft = 1000;
end

NHarm = size(X{1},1)-1;

for i = 1:length(X)
    x{i} = permute(freq2time(X{i},NHarm,Nfft),[2 3 1]);
    x{i} = cat(3,x{i},x{i}(:,:,1));
    amp{i} = max(x{i},[],3) - min(x{i},[],3);
end

if length(x) == 1
    x = x{1};
    amp = amp{1};
end