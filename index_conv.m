function x = index_conv(k, n)
k = k-1;
x(2) = mod(k, n)+1;
x(1) = floor(k/n)+1;
