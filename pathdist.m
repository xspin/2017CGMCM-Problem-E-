function d = pathdist(TABC, path, k)
d = 0;
for i=1:length(path)-1
    d = d + TABC{k}(path(i), path(i+1));
end