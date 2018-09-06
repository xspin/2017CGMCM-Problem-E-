function path = spath(dist, trace, ij)
% return shortest path
i = ij(1);
j = ij(2);
path = [i];
while i ~= j
    path = [path trace(i, j)];
    i = trace(i,j);
end
