function [Ared, Cred, dAred, dCred] = hbm_reduce(hbm,problem,A,dA)
if problem.NNL == problem.NDof
    Ared = A;
    Cred = speye(hbm.harm.NRetain);
    if nargout > 2
        dAred = dA;
        dCred = spalloc(hbm.harm.NRetain,hbm.harm.NRetain,0);
    end
else
    Cred  = zeros(hbm.harm.NRetainNL,hbm.harm.NRetain);
    dCred = zeros(hbm.harm.NRetainNL,hbm.harm.NRetain);

    iLin = hbm.harm.iLin;
    iNL = hbm.harm.iNL;

    App = A(iLin,iLin);
    Aqp = A(iNL,iLin);
    Apq = A(iLin,iNL);
    Aqq = A(iNL,iNL);
    App_inv = inv(App);
    
    Ared(iNL,:) = Aqq;
    Ared(iLin,:)= Apq;
    
    Cred(:,iLin) = -Aqp*App_inv;
    Cred(:,iNL) = eye(length(iNL));

    if nargout > 2
        dApp = dA(iLin,iLin);
        dAqp = dA(iNL,iLin);
        dApq = dA(iLin,iNL);
        dAqq = dA(iNL,iNL);
        dApp_inv = -App_inv*dApp*App_inv;
        
        dAred(iNL,:)  = dAqq;
        dAred(iLin,:) = dApq;
    
        dCred(:,iLin) = -dAqp*App_inv-Aqp*dApp_inv;
        dCred(:,iNL)  = 0;
    end
end