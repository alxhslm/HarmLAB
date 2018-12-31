function varargout = hbm_constraints(part,hbm,problem,w0,x,u)

if ~iscell(part)
    part = {part};
end

for i = 1:length(part)
    switch part{i}
        case 'func'
            varargout{i} = w0(1) - x(3);
        case 'jacobX'
            Jx = 0*x';
            Jx(3) = -1;
            varargout{i} = Jx;
        case 'derivW'
            Dw = 1;
            varargout{i} = Dw;
        case 'derivA'
            Da = 0;
            varargout{i} = Da;
    end
end