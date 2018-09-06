function [DFZF, ds, Fe] = select(traceZ, tDF1, I1, tDF2, I2, n, Fe)
% DFZF: selected paths  
% ds: distances
% Fs: selected F nodes
global index_J index_D index_Z index_F n_Z;
DFZF = [];
ds = [];
%Fe = [];
n_F = size(traceZ, 1);

i1 = (find(tDF1, 1)); % index of the 1st non-zero element
i2 = (find(tDF2, 1));

% limit = ones(1,n_Z)*sum(n)/n_Z+1;
while any(n)
    p1 = index_conv(I1(i1), n_F);
    while ~(isempty(find(p1(1)==Fe,1)) && isempty(find(p1(2)==Fe,1)))
        i1 = i1+1;
        p1 = index_conv(I1(i1), n_F);
    end
    p2 = index_conv(I2(i2), n_F);
    while ~(isempty(find(p2(1)==Fe,1)) && isempty(find(p2(2)==Fe,1)))
        i2 = i2+1;
        p2 = index_conv(I2(i2), n_F);
    end
    z = find(traceZ(p1(1), p1(2))==index_Z);
%     if limit(z) <= 0
%         if all(n>0)
%             if tDF1(i1)<tDF2(i2)
%                 i1 = i1+1;
%             else
%                 i2=i2+1;
%             end
%         elseif n(2)==0
%             i1 = i1+1;
%         elseif n(1)==0
%             i2=i2+1;
%         end
%         continue;
%     end
%     limit(z) = limit(z)-1
    
    
    if all(n>0)
        if tDF1(i1)<tDF2(i2)
            ds = [ds tDF1(i1)];
%             DFZF = [DFZF; idx('D01'), idx('F01')-1+p1(1), traceZ(p1(1), p1(2)), idx('F01')-1+p1(2)];
            DFZF = [DFZF; index_D(1), index_F(p1(1)), traceZ(p1(1), p1(2)), index_F(p1(2))];
            Fe = [Fe, p1(1), p1(2)];
            i1 = i1+1;
            n(1) = n(1)-1;
        else
            ds = [ds tDF2(i2)];
%             DFZF = [DFZF; idx('D02'), idx('F01')-1+p2(1), traceZ(p2(1), p2(2)), idx('F01')-1+p2(2)];
            DFZF = [DFZF; index_D(2), index_F(p2(1)), traceZ(p1(1), p1(2)), index_F(p2(2))];
            Fe = [Fe, p2(1), p2(2)];
            i2 = i2+1;
            n(2) = n(2)-1;
        end
    elseif n(2)==0
        ds = [ds tDF1(i1)];
%         DFZF = [DFZF; idx('D01'), idx('F01')-1+p1(1), traceZ(p1(1), p1(2)), idx('F01')-1+p1(2)];
          DFZF = [DFZF; index_D(1), index_F(p1(1)), traceZ(p1(1), p1(2)), index_F(p1(2))];
        Fe = [Fe, p1(1), p1(2)];
        i1 = i1+1;
        n(1) = n(1)-1;
    elseif n(1)==0
        ds = [ds tDF2(i2)];
%         DFZF = [DFZF; idx('D02'), idx('F01')-1+p2(1), traceZ(p2(1), p2(2)), idx('F01')-1+p2(2)];
          DFZF = [DFZF; index_D(2), index_F(p2(1)), traceZ(p1(1), p1(2)), index_F(p2(2))];
        Fe = [Fe, p2(1), p2(2)];
        i2 = i2+1;
        n(2) = n(2)-1;
    end
end