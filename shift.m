function p = shift(p, d, tag)
if tag == 1
    return
end
for i=1:size(p,1)-1
    dx = p(i,1)-p(i+1,1);
    dy = p(i,2)-p(i+1,2);
    theta = pi/2-atan(dy/dx);
    p(i,1) = p(i,1) + cos(theta);
    p(i,2) = p(i,2) + sin(theta);
end