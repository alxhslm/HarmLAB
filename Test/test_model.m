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
             fel = P.knl .* sgn_power(P.R*x,P.n);
             fdamp = P.cnl .* (P.R*xdot).^3;
             Fnl = P.R' * (fel + fdamp);
             
             varargout{end+1} = Fnl;
        case 'nl_x'   
             [~,d] = sgn_power(P.R*x,P.n);
             Knl = zeros(problem.NDof);
             for j = 1:size(P.R,1)
                Knl = Knl + mtimesx(P.knl(j) .* permute(d(j,:),[1 3 2]),P.R(j,:)'*P.R(j,:));
             end
             varargout{end+1} = Knl;
         case 'nl_xdot'   
             Cnl = zeros(problem.NDof);
             d = 3 * P.cnl .* (P.R*xdot).^2;
             for j = 1:size(P.R,1)
                Cnl = Cnl + mtimesx(P.cnl(j) .* permute(d(j,:),[1 3 2]),P.R(j,:)'*P.R(j,:));
             end
             varargout{end+1} = Cnl;
        case 'nl_u'
            varargout{end+1} = zeros(2,1,NPts);
        case 'nl_udot'
            varargout{end+1} = zeros(2,1,NPts);
        otherwise
            varargout{end+1} = [];
    end
end