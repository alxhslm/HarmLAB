function [Bk,Bc,Bm,Bx] = linear_jacobian(harm,Ku,Cu,Mu)

%terms dependent on each frequency
Ak{1} = Ku;
for j = 2:harm.NFreq
    Ak{j} = blkdiag(Ku,Ku);
end

for k = 1:2
    Ac{k}{1} = 0*Cu;
    Am{k}{1} = 0*Mu;
    for j = 2:harm.NFreq
        Ac{k}{j} =  harm.rFreqBase(k)*harm.kHarm(j,k)*antiblkdiag(-Cu,Cu);
        Am{k}{j} = -(harm.rFreqBase(k)*harm.kHarm(j,k))^2*blkdiag(Mu,Mu);
    end
end

Bk = blkdiag(Ak{:});
for k = 1:2
    Bc{k} = blkdiag(Ac{k}{:});
    Bm{k} = blkdiag(Am{k}{:});
end

%cross-terms
Bx{1} = 0*Mu;
for j = 2:harm.NFreq
    Bx{j} = -2*prod(harm.kHarm(j,:).*harm.rFreqBase)*blkdiag(Mu,Mu);
end
Bx = blkdiag(Bx{:});

Bk = sparse(Bk);
for k = 1:2
    Bm{k} = sparse(Bm{k});
    Bc{k} = sparse(Bc{k});
end
Bx = sparse(Bx);
