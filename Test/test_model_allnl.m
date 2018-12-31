function varargout = test_model(part,t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0)

P = problem.P;

switch part
    case 'excite'
        U = zeros(hbm.harm.NFreq,2);
        U(1,2) = P.f0;
        U(2,2) = P.f;
        if hbm.harm.NHarm(2)>0
            ii = hbm.harm.kHarm(:,1) == 0 & hbm.harm.kHarm(:,2) == 1; 
            U(ii,2) = P.f2;
        end
        varargout= {U};
    case  {'nonlin','output'}
         Fe = problem.Ku*u + problem.Cu*udot + problem.Mu*uddot;
         Fl = -problem.C*xdot - problem.K*x - problem.M*xddot;
         x1 = x(1,:); x2 = x(2,:);
         f = sgn_power(x1-x2,P.n);
         Fnl = [-1;1]*P.knl * f - P.cnl * xdot.^3;
         varargout{1} = Fl + Fnl + Fe;
    case  'statederiv'        
         x1 = permute(x(1,:),[1 3 2]);
         x2 = permute(x(2,:),[1 3 2]);
         xdot1 = permute(xdot(1,:),[1 3 2]);
         xdot2 = permute(xdot(2,:),[1 3 2]);
         [~,d] = sgn_power(x1-x2,P.n);
         k = P.knl * d;
         c1 = P.cnl * 3 * xdot1.^2;
         c2 = P.cnl * 3 * xdot2.^2;
         
         Knl = -[k -k; -k k];
         Cnl = -[c1,0*c1;0*c2,c2];
                  
         varargout = {[],[]};
    case  'inputderiv'
        varargout = {zeros(2,2),zeros(2,2)};
    case 'freqderiv'
        varargout{1} = [];
end