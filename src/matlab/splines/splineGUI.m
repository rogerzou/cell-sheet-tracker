function splineGUI(spl, img)

if nargin < 2 || isempty(img)
    state.image = false;
else
    [x y] = round(size(img));
    state.image = true;
end

state.d = spl;
state.np = size(state.d.control, 2);

if state.image
    state.scale = 20;
    state.d = splineAffineXform(state.d, (state.scale - 1) * eye(2), [x; y]);
    
    if ~exist('grad.m', 'file')
        upath = sprintf('%cUsers%ccarlo%cDropbox%cDocuments%cCode%cmatlab%cUTILS', ...
            filesep, filesep, filesep, filesep, filesep, filesep, filesep);
        addpath(genpath(upath))
    end
    
    % Compute the magnitude of the image gradient
    sigma = 2;
    state.gradImg = grad(img, sigma);
    state.gradImg = sqrt(state.gradImg{1}.^2 + state.gradImg{2}.^2);
else
    state.scale = 1;
end

hold off
if state.image
    imshow(img);
    hold on
end
state.handles = splineDraw(state.d);
axis ij
axis off
state = redraw(state);

state.pointIdx = 0;
state.cp = NaN;

zoom off;

state.WindowButtonDownFcn = get(gcf, 'WindowButtonDownFcn');
state.WindowButtonMotionFcn = get(gcf, 'WindowButtonMotionFcn');
state.WindowButtonUpFcn = get(gcf, 'WindowButtonUpFcn');
state.CloseRequestFcn = get(gcf, 'CloseRequestFcn');
state.KeyPressFcn = get(gcf, 'KeyPressFcn');

set(gcf, 'WindowButtonDownFcn', @selectPoint);
set(gcf, 'WindowButtonMotionFcn', @trackPoint);
set(gcf, 'WindowButtonUpFcn', @stopTracking);
set(gcf, 'CloseRequestFcn', @closeFigure);
set(gcf, 'KeyPressFcn', @cleanup);

% Save the state with the current axes, so different axes can have
% different GUIs
set(gca, 'UserData', state);

    function selectPoint(~, ~)
        s = get(gca, 'UserData');
        s.cp = get(gca,'CurrentPoint');
        s.cp = s.cp(1, 1:2)';
        dist = s.d.control - repmat(s.cp, [1 s.np]);
        dist = sum(dist .^ 2, 1);
        [~, s.pointIdx] = min(dist);
        set(gca, 'UserData', s);
    end

    function trackPoint(~, ~)
        s = get(gca, 'UserData');
        if s.pointIdx
            newcp = get(gca,'CurrentPoint');
            newcp = newcp(1, 1:2)';
            if norm(newcp - s.cp) > 0.001 * s.scale
                s.cp = newcp;
                s.d.control(:, s.pointIdx) = s.cp;
                if s.d.open && (s.pointIdx == 1 || s.pointIdx == s.np)
                    % First and last point do not move
                    return
                end
                s.d = splineEvalEven(s.d);
                s = redraw(s);
                set(gca, 'UserData', s);
            end
        end
    end

    function stopTracking(~, ~)
        s = get(gca, 'UserData');
        s.pointIdx = 0;
        set(gca, 'UserData', s);
    end

    function s = redraw(s)
        splineDraw(s.d, s.handles);
        if s.image
            title(matchQuality(s.d, s.gradImg))
        end
    end

    function cleanup(~, ~)
        s = get(gca, 'UserData');
        set(gca, 'UserData', s);
        
%         % Save splines from all subplots of this figure
%         f = gcf;
%         child = get(get(gca, 'Parent'), 'Children');
%         for c = 1:length(child)
%             h = child(c);
%             s = get(h, 'UserData');
%             if ~isempty(s) && isstruct(s) && isfield(s, 'd')
%                 filename = sprintf('savedSpline.%d.%d.mat', f, c);
%                 d = s.d; %#ok<NASGU>
%                 save(filename, 'd');
%                 fprintf(1, 'Saved spline from figure %d, subplot %d in file %s\n', ...
%                     f, c, filename);
%             end
%         end

        zoom off;
        
        set(gcf, 'WindowButtonDownFcn', s.WindowButtonDownFcn);
        set(gcf, 'WindowButtonMotionFcn', s.WindowButtonMotionFcn);
        set(gcf, 'WindowButtonUpFcn', s.WindowButtonUpFcn);
        set(gcf, 'CloseRequestFcn', s.CloseRequestFcn);
        set(gcf, 'KeyPressFcn', s.KeyPressFcn);
    end

    function closeFigure(~, ~)
        cleanup;
        delete(gcf);
    end
end