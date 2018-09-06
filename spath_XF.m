function [XF, Fe] = spath_XF(Dist, ks, n, Fe)
%
global index_F;
XF = [];
if n<1
    return 
end
d = [];
for k=ks
    for i=index_F
        if isempty(find(i==Fe,1))
            d = [d, Dist(k, i)];
            XF = [XF; k i];
        end
    end
end
[~, I] = sort(d);
XF = XF(I(1:n),:);
Fe = [Fe, XF(:,2)'];