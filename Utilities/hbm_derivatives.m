function varargout = hbm_derivatives(fun,var,States,hbm,problem)
NPts = size(States.t,2);

if ~iscell(var)
    var = {var};
end

for i = 1:length(var)
    if strcmp(var{i},'w')    
        %for frequency derivatives, it's a bit more complicated as we have to rescale the time vector as well     
        df_dw = feval(problem.model,[fun '_w'],States,hbm,problem);
        if isempty(df_dw)
            h = 1E-6;
            NFreq = length(States.w0);
            
            if hbm.dependence.w
                for k = 1:NFreq
                    States2 = States;
                    States2.w(k) = States.w(k) + h;
                    States2.t(k,:) = w0(k)*States2.t(k,:)/States2.w(k);
                    States2.f = feval(problem.model,fun,States2,hbm,problem);
                    df_dw{k} = (States2.f-States.f)./h;
                end
            else
                df_dw = repmat({zeros(size(f,2),NPts)},1,NFreq);
            end
        end
        for k = 1:length(w0)
            df_dw{k} = df_dw{k}.';
        end
        varargout{i} = df_dw;
    else
        %for everything else, we just peturb each element in turn
        df_dvar = feval(problem.model,[fun '_' var{i}],States,hbm,problem);
        if isempty(df_dvar)
            if hbm.dependence.x
                h = 1E-10;
                df_dvar = zeros(size(States.f,1),size(States.(var{i}),1),NPts);
                for j = 1:problem.NDof
                    States2 = States;
                    States2.(var{i})(j,:) = States.(var{i})(j,:) + h;
                    States2.f = feval(problem.model,fun,States2,hbm,problem);
                    df_dvar(:,j,:) = permute((States2.f-States.f)./h,[1 3 2]);
                end
            else
                df_dvar = zeros(size(States.f,2),size(States.(var{i}),1),NPts);
            end
        end
        varargout{i} = df_dvar;
    end
end