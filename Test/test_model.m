function varargout = test_model(part,States,hbm,problem)
if ~iscell(part)
    part = {part};
end

P = problem.P;
NPts = size(States.t,2);
x = States.x;
xdot = States.xdot;

varargout = {};

for i = 1:length(part)
    switch part{i}           
        %% Nonlin
        case  {'nl','output'}
             x1 = x(1,:); x2 = x(2,:);
             f = sgn_power(x1-x2,P.n);
             Fnl = [1;-1]*P.knl * f + P.cnl * xdot(1:2,:).^3;
             Fjen = 0;
             
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
             
             varargout{end+1} = Knl;
         case 'nl_xdot'   
             xdot1 = permute(xdot(1,:),[1 3 2]);
             xdot2 = permute(xdot(2,:),[1 3 2]);
             c1 = P.cnl * 3 * xdot1.^2;
             c2 = P.cnl * 3 * xdot2.^2;

             Cnl = [c1,0*c1;0*c2,c2];
             varargout{end+1} = Cnl;
        case 'nl_u'
            varargout{end+1} = zeros(2,1,NPts);
        case 'nl_udot'
            varargout{end+1} = zeros(2,1,NPts);
        otherwise
            varargout{end+1} = [];
    end
end