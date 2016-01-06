% mod operation consistent with Matlab style indexing (starting at 1
% instead of 0

function y = mmod(x, m)

y = mod(x - 1, m) + 1;