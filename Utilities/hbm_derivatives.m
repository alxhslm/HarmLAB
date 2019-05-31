function varargout = hbm_derivatives(fun,var,States,hbm,problem)
NPts = size(States.t,2);

if ~iscell(var)
    var = {var};
end

for i = 1:length(var)
    switch var{i}(1)
        case 'w'
            %for frequency derivatives, it's a bit more complicated as we have to rescale the time vector as well
            df_dw = feval(problem.model,[fun '_w'],States,hbm,problem);
            if isempty(df_dw)
                h = 1E-6;
                NFreq = length(States.w0);
                
                if hbm.dependence.w
                    for k = 1:NFreq
                        States2 = States;
                        States2.w0(k) = States.w0(k) + h;
                        States2.t(k,:) = States.w0(k)*States2.t(k,:)/States2.w0(k);
                        States2.f = feval(problem.model,fun,States2,hbm,problem);
                        df_dw{k} = (States2.f-States.f)./h;
                    end
                else
                    df_dw = repmat({zeros(size(f,1),NPts)},1,NFreq);
                end
            end
            for k = 1:length(States.w0)
                df_dw{k} = df_dw{k}.';
            end
            varargout{i} = df_dw;
        case 'x'
            %for everything else, we just peturb each element in turn
            df_dx = feval(problem.model,[fun '_' var{i}],States,hbm,problem);
            if isempty(df_dx)
                if hbm.dependence.(var{i})
                    h = 1E-10;
                    df_dx = zeros(size(States.f,1),problem.NDof,NPts);
                    for j = problem.iNL'
                        States2 = States;
                        States2.(var{i})(j,:) = States.(var{i})(j,:) + h;
                        States2.f = feval(problem.model,fun,States2,hbm,problem);
                        df_dx(:,j,:) = permute((States2.f-States.f)./h,[1 3 2]);
                    end
                else
                    df_dx = zeros(size(States.f,1),problem.NDof,NPts);
                end
            end
            varargout{i} = df_dx;
        case 'u'
            %for everything else, we just peturb each element in turn
            df_du = feval(problem.model,[fun '_' var{i}],States,hbm,problem);
            if isempty(df_du)
                if hbm.dependence.(var{i})
                    h = 1E-10;
                    df_du = zeros(size(States.f,1),problem.NInput,NPts);
                    for j = 1:problem.NInput
                        States2 = States;
                        States2.(var{i})(j,:) = States.(var{i})(j,:) + h;
                        States2.f = feval(problem.model,fun,States2,hbm,problem);
                        df_du(:,j,:) = permute((States2.f-States.f)./h,[1 3 2]);
                    end
                else
                    df_du = zeros(size(States.f,1),problem.NInput,NPts);
                end
            end
            varargout{i} = df_du;
        otherwise
            error('Unrecognized input')
    end
end