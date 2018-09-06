function v = variance(x)
n = size(x,1)-1;
x = x-mean(x);
x = (sum(x.^2));
v = sqrt(x/n);
v = sum(v);