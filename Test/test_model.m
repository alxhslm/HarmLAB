function varargout = test_model(part,t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0)
if ~iscell(part)
    part = {part};
end

P = problem.P;
NPts = size(x,2);

varargout = {};

for i = 1:length(part)
    switch part{i}           
        %% Nonlin
        case  {'nl','output'}
             x1 = x(1,:); x2 = x(2,:);
             f = sgn_power(x1-x2,P.n);
             Fnl = [1;-1]*P.knl * f + P.cnl * xdot(1:2,:).^3;

             if problem.NAlg>0
                 f = P.kJen*(x1-x2+xalg);
                 Fjen = [1;-1]*f;
             else
                 Fjen = 0;
             end
             varargout{end+1} = Fnl + Fjen;
        case 'nl_x'   
            varargout{end+1} = [];
            return
             x1 = permute(x(1,:),[1 3 2]);
             x2 = permute(x(2,:),[1 3 2]);
             [~,d] = sgn_power(x1-x2,P.n);
             k = P.knl * d;

             Knl = [k -k;
                   -k  k];
             
             if problem.NAlg>0
                 Kalg = P.kJen*[1 -1;
                               -1  1];
                 Kalg = repmat(Kalg,1,1,NPts);
             else
                 Kalg = 0;
             end
             
             varargout{end+1} = Knl + Kalg;
         case 'nl_xdot'   
             xdot1 = permute(xdot(1,:),[1 3 2]);
             xdot2 = permute(xdot(2,:),[1 3 2]);
             c1 = P.cnl * 3 * xdot1.^2;
             c2 = P.cnl * 3 * xdot2.^2;

             Cnl = [c1,0*c1;0*c2,c2];
             varargout{end+1} = Cnl;
        case 'nl_u'
            varargout{end+1} = zeros(2,1,length(t));
        case 'nl_udot'
            varargout{end+1} = zeros(2,1,length(t));
         case 'nl_alg'   
             if problem.NDof > 0
                 Kalg = P.kJen*[1; -1];
                 Kalg = repmat(Kalg,1,1,NPts);
                 varargout{end+1} = Kalg;
             else
                 varargout{end+1} = [];
             end
            
        %% Algebraic
        case 'alg'
            if problem.NDof > 0
                x1 = x(1,:); x2 = x(2,:);
                xdot1 = xdot(1,:); xdot2 = xdot(2,:);
                f = P.kJen*(x1-x2+xalg);
                varargout{end+1} = f - P.FJen*sgnSmooth(xdot1-xdot2,1E-4);
            else
                varargout{end+1} = [];
            end
        otherwise
            varargout{end+1} = [];
    end
end