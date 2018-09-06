function [distance, trace] = Floyd(edges)
% Shortest path problem
% input: distance matrix
% output: distance and shortest path

trace = [];
distance = edges;
n = size(distance,1);

for i = 1:n
    for j = 1:n
        trace(i,j) = j;   
    end
end

for k = 1:n
    for i = 1:n
        for j = 1:n
            if distance(i,k) + distance(k,j) < distance(i,j)
                distance(i,j) = distance(i,k) + distance(k,j);
                trace(i,j) = trace(i,k);
            end
        end
    end
end
