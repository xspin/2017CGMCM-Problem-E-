function k = vidx(t, i)
% n = [6 6 12];
global n_ABC;
n = sum(n_ABC');
if nargin==2
    k = sum(n(1:t-1))+i;
else
    k = 1 + (t>n(1)) + (t>n(1)+n(2));
end