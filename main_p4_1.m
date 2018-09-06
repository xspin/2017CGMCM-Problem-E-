% 问题四 主程序 计算多费程与毁坏发射点数目
% 打击单个道路节点
clc;
clear;

fid = load('A.mat');
A = single(cell2mat(struct2cell(fid)));

[r,c] = size(A);
delta_S = zeros(62,1); % 多费程
num_F = zeros(62,1);   % 破坏发射点数

h = waitbar(0);
for n = 1:62           % J01 ~ J62
    %炸毁n点后新的邻接矩阵
    waitbar(n/62,h,['Now processing ' num2str(n) ' / ' num2str(62)]);
    
    A_new = A;
    A_new(n, :) = inf;
    A_new(:, n) = inf;
    A_new(n, n) = 0;
    
    ind = find(A(n,:) ~= inf & A(n,:) ~= 0);
    ind2 = ind(ind <= 70);             % 找到与n相邻的节点的编号
    num_F(n) = length(ind(ind >= 71)); % 大于等于71的为发射点数编号
    
    if length(ind2) >= 2
        comb = nchoosek(ind2,2);       % 下标两两组合
        for i = 1:size(comb,1)
            [dist1] = Floyd(A,comb(i,1),comb(i,2));      % 估计未炸毁前两点间的最短距离和路径
            [dist2] = Floyd(A_new,comb(i,1),comb(i,2));  % 估计炸毁后两点间新的最短距离和路径
            delta_S(n) = delta_S(n) + (dist2 - dist1);   % 计算多费程
        end
    end
end
close(h)

S_all = sum(A(find(A ~= inf & A ~= 0))) ./ 2;
lambda = (num_F .* ( num_F + 1) .* (2 * num_F + 1) ./ 6).^ 0.5 ./...
    (60 * (60 + 1) * (2 * 60 + 1) / 6) .^ 0.5;

S_k = S_all * lambda;
Eff = S_k + delta_S;

[Eff, index] = sort(Eff, 'descend');
delta_S = delta_S(index);
num_F = num_F(index);