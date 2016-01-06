function [q, qgrad] = matchQuality(d, g)

q = quality(d.curve);

if nargout == 2
    delta = 1;
    qgrad = zeros(size(d));
    for i = 1:size(d.control, 1)
        for j = 1:size(d.control, 2)
            d1 = d;
            d1.control(i, j) = d1.control(i, j) + delta;
            d1 = splineEval(d1);
            q1 = quality(d1.curve);
            qgrad(i, j) = (q1 - q) / delta;
        end
    end
end

    function qq = quality(bb)
        % No harm in closing a closed curve, since deltaT is zero for the
        % last interval in that case
        wbb = [bb bb(:, 1)];
        deltaT = sqrt(sum(diff(wbb, 1, 2) .^ 2, 1));
        qq = sum(interp2(g, bb(1, :), bb(2, :)) .* deltaT);
    end

end