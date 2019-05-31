function [Bk,Bc,Bm,Bx] = linear_jacobian(harm,Ku,Cu,Mu)

%terms dependent on each frequency
Ak{1}{1} = Ku;
Ac{1}{1} = 0*Cu;
Am{1}{1} = 0*Mu;

Ak{2}{1} = 0*Ku;
Ac{2}{1} = 0*Cu;
Am{2}{1} = 0*Mu;

for j = 2:harm.NFreq
    Ak{1}{j} = blkdiag(Ku,Ku);
    Ak{2}{j} = blkdiag(0*Ku,0*Ku);
end

for k = 1:2
    for j = 2:harm.NFreq
        Ac{k}{j} =  harm.rFreqBase(k)*harm.kHarm(j,k)*antiblkdiag(-Cu,Cu);
        Am{k}{j} = -(harm.rFreqBase(k)*harm.kHarm(j,k))^2*blkdiag(Mu,Mu);
    end
end

for k = 1:2
    Bk{k} = blkdiag(Ak{k}{:});
    Bc{k} = blkdiag(Ac{k}{:});
    Bm{k} = blkdiag(Am{k}{:});
end

%cross-terms
Bx{1} = 0*Mu;
for j = 2:harm.NFreq
    Bx{j} = -2*prod(harm.kHarm(j,:).*harm.rFreqBase)*blkdiag(Mu,Mu);
end
Bx = blkdiag(Bx{:});

for k = 1:2
    Bk{k} = sparse(Bk{k});
    Bm{k} = sparse(Bm{k});
    Bc{k} = sparse(Bc{k});
end
Bx = sparse(Bx);
