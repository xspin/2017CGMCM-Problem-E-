function h = drawArrow(fid, xs, ys, varargin)

for i=1:length(xs)-1
    x = xs(i:i+1);
    y = ys(i:i+1);
    h = quiver(fid, x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,varargin{:});
end
