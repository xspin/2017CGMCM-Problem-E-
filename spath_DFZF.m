function [DF1, DF2, traceZ] = spath_DFZF(Dist, alpha, beta, gamma)
% DFi: shortest distance/path from Di,..., Fj,..., Zk ,... to Fl.
% traceZ: traceZ(i,j) is the node index of Zt in the path DFk(i,j)

global n_J  n_D n_Z n_F;
global index_J index_D index_Z index_F;

DF1 = [];   % shortest paths starting from D1
DF2 = [];   % shortest paths starting from D2
%DF = [];    % shortest paths without the D1/D2 limitation
% DFij ~= DFji
% n_F = 60;
%traceD = zeros(n_F, n_F, 'uint32');
traceZ = zeros(n_F, n_F, 'uint32');% all are same
% alpha = 0.0005;
% beta = 0.1;
% gamma = 0.1;
for i=1:n_F
    for j=1:n_F
        if i==j
            continue
        end
%         i_D = idx('D01')-1+[1:n_D];
%         i_Z = idx('Z01')-1+[1:n_Z];
%         i_F1 = idx('F01')-1+i;
%         i_F2 = idx('F01')-1+j;
        i_F1 = index_F(i);
        i_F2 = index_F(j);
        d1 = Dist(i_F1, index_D)*alpha;    
        [d2, tz] = min(Dist(i_F1, index_Z)*beta+Dist(i_F2, index_Z)*gamma);
        DF1(i,j) = d1(1)+d2;
        DF2(i,j) = d1(2)+d2;
        %DF(i,j) = min(DF1(i,j), DF2(i,j));
        %td = (d1(1)<d1(2))+1;
        %traceD(i,j) = i_D(td);
        traceZ(i,j) = index_Z(tz);
    end
end