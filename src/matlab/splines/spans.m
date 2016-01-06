function span = spans(knot, order, open)

span = find(knot(1:(end-1)) < knot(2:end));
span = span(span >= order & span <= length(knot) - order);

if ~open
    span = [span length(knot) - order + 1];
end