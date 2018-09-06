% 问题四 主程序 计算多费程与毁坏发射点数目
% 以同时打击两个道路节点为例
clc;
clear;

fid = load('A.mat');
A = single(cell2mat(struct2cell(fid)));

[S_k, delta_S, num_F ] = skDeltas( fid );

comb = nchoosek(1:62,2);
temp = find(A ~= inf & A ~= 0);
S_all = sum(A(temp)) ./ 2;

lambda_pair = zeros(length(comb(:,1)),1);
num_F_pair = zeros(length(comb(:,1)),1);
delta_S_pair = zeros(length(comb(:,1)),1);
S_k_pair = zeros(length(comb(:,1)),1);

h = waitbar(0);
for i = 1:length(comb(:,1))
    waitbar(i/length(comb(:,1)),h,['Now processing ' num2str(i) ' / ' num2str(length(comb(:,1)))]);
    
    if comb(i,1) ~= 2 && comb(i,2) ~= 2
        delta_S_pair(i) = sum(delta_S(comb(i, :)));
        num_F_pair(i) = sum(num_F(comb(i, :)));
        lambda_pair(i) = (num_F_pair(i) .* (  num_F_pair(i) + 1) .* (2 *  num_F_pair(i)+ 1) ./ 6).^ 0.5 ./...
            (60 * (60 + 1) * (2 * 60 + 1) / 6) .^ 0.5;
        S_k_pair(i) = S_all .* lambda_pair(i);
        
        %%% if using two-dimensional parameters,please add 'm(i)' and 'n(i)'
        
        m(i) = comb(i, 1);
        n(i) = comb(i, 2);

    end
end
close(h)

dbstop if error

Eff_pair = S_k_pair + delta_S_pair;
[Eff_pair, index] = sort(Eff_pair, 'descend');
delta_S_pair = delta_S_pair(index);
num_F_pair = num_F_pair(index);
index = [m(index)',n(index)'];