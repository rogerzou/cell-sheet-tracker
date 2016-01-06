function h = splineDraw(s, h, drawNeedles)

if nargin < 3 || isempty(drawNeedles)
    drawNeedles = false;
end

control = s.control;
if ~s.open
    control = [control control(:, 1)];
end

if drawNeedles
    nx = [s.needle.in(1, :); s.needle.out(1, :)];
    ny = [s.needle.in(2, :); s.needle.out(2, :)];
end

if nargin == 2
    set(h.curve, 'XData', s.curve(1, :), 'YData', s.curve(2, :));
    set(h.ctrlLine, 'XData', control(1, :), 'YData', control(2, :));
    set(h.ctrlPoint, 'XData', control(1, :), 'YData', control(2, :));
    if drawNeedles
        for n = 1:length(h.needle)
            set(h.needle(n), 'XData', nx(:, n), 'YData', ny(:, n));
        end
    end
else
    oldhold = ishold;
    h.curve = plot(s.curve(1, :), s.curve(2, :), '-r', 'LineWidth', 4);
    hold on
    h.ctrlLine = plot(control(1, :), control(2, :), '-y');
    h.ctrlPoint = plot(control(1, :), control(2, :), '.y', 'MarkerSize', 12);
    if drawNeedles
       h.needle = plot(nx, ny, '-g');
    end
    axis equal
    if ~oldhold
        hold off
    end
    axis ij
end