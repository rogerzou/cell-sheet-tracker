function f = displayGraph(I, V, E, visible)
%DISPLAYGRAPH Displays graph with splines and vertices
% 
% INPUTS
% I: (2D matrix) image.
% V: cell array of 2D vectors of vertices.
% Splines: cell array of spline structs.
% visible: {'on (default)', 'off'} whether the figure is visible.
% 
% OUTPUTS
% f: figure handle.
% 
% @author Roger Zou
% @date 5/19/15

% validate input
if nargin < 4 || isempty(visible)
    visible = 'on';
end

% get lengths
N = length(V);
M = length(E);

% setup figure
f = figure('Visible', visible);
hold on;
imagesc(I);
colormap gray

% draw the vertices
for ii=1:N

    % get current vertex
    v_i = V{ii};

    % draw the vertices
    if ~sum(isnan(v_i))
        plot(v_i(1), v_i(2),'go',...
            'MarkerSize',5,...
            'MarkerEdgeColor','y',...
            'MarkerFaceColor','y');
    end
    
end

% draw the splines
for ii=1:M

    % get current spline
    si = E{ii};

    % draw the spline and the sample points on it
    if ~isempty(si)
        line(si.curve(1,:), si.curve(2,:), 'Color', 'y', 'LineWidth', 2)
%         line(si.control(1, :), si.control(2, :), ...
%             'Marker', '.', 'MarkerSize', 10, 'Color', 'y');
    end
    
end

axis equal;
axis xy;

end
