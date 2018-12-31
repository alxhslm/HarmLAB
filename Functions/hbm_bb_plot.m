function varargout = hbm_bb_plot(command,hbm,problem,a,h,w)
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
            
            A = a;
            W = w;
            H = h;
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
                A(end+1) = a;
                H(end+1) = h;
                W(end+1) = w;
                [A,isort] = sort(A);
                H = H(isort); W = W(isort);
                set(hBB(1),'xdata',A,'ydata',W);
                set(hBB(2),'xdata',A,'ydata',H);

                %update our progress
                if ~ishandle(hWaitbar)
                     hWaitbar = waitbar(0, 'Amplitude Range');
                end
                waitbar((a-problem.A0)/(problem.AEnd - problem.A0),hWaitbar);
            else
                plot(axBB(1),a,w,'r.')
                plot(axBB(2),a,h,'r.')
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
matlabPos = getMatlabSize;
figPos = matlabPos;
figPos(4) = matlabPos(4)/2;
figPos(2) = matlabPos(2) + figPos(4);
fig = figure('Name',[problem.name],'OuterPosition',figPos,'WindowStyle', 'Docked');

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