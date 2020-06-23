function varargout = hbm_bb_plot(command,hbm,problem,results)
persistent figBB axBB hBB hWaitbar A H W
if hbm.cont.bUpdate
    switch command
        case 'init'
            if ~isempty(figBB) && ishandle(figBB)
                close(figBB)
                hBB = [];
            end
            if ~isempty(hWaitbar) && ishandle(hWaitbar)
                close(hWaitbar)
                hWaitbar = [];
            end
            
            A = results.A;
            W = results.w;
            H = abs(results.H);
            [figBB,axBB,hBB] = createFig(hbm,problem,A,H,W);        
            hWaitbar = waitbar(0, 'Amplitude Range');
            if nargout > 0
                varargout{1} = figBB;
                if  nargout > 1
                    varargout{2} = axBB;
                end
            end

        case {'data','err'}
            if ~ishandle(figBB)
                [figBB,axBB,hBB] = createFig(hbm,problem,A,H,W);
            end
            if strcmp(command,'data')
                A(end+1) = results.A;
                H(end+1) = abs(results.H);
                W(end+1) = results.w;
                [A,isort] = sort(A);
                H = H(isort); W = W(isort);
                set(hBB(1),'xdata',A,'ydata',W);
                set(hBB(2),'xdata',A,'ydata',H);

                %update our progress
                if ~ishandle(hWaitbar)
                     hWaitbar = waitbar(0, 'Amplitude Range');
                end
                waitbar((results.A-problem.A0)/(problem.AEnd - problem.A0),hWaitbar);
            else
                plot(axBB(1),results.A,results.w,'r.')
                plot(axBB(2),results.A,abs(results.H),'r.')
            end
            drawnow
        case 'close'
            close(hWaitbar)
            close(figBB)
            hWaitbar = []; 
            figBB = [];         
            hBB = [];
            axBB = [];
    end
end


function [fig,axBB,hBB] = createFig(hbm,problem,A,H,w)
fig = figure('Name',['BB: ' problem.name]);

axBB(1) = subplot(2,1,1);
hBB(1) = plot(axBB(1),A,w,'g.-');
hold on
ylabel('Frequency (rad/s)')
xlabel('Amplitude')
xlim(axBB(1),[problem.AMin problem.AMax]);

axBB(2) = subplot(2,1,2);
hBB(2) = plot(axBB(2),A,H,'g.-');
hold on
ylabel('Peak amp')
xlabel('Amplitude')
xlim(axBB(2),[problem.AMin problem.AMax]);