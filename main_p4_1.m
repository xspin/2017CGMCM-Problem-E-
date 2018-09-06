% ������ ������ �����ѳ���ٻ��������Ŀ
% ���������·�ڵ�
clc;
clear;

fid = load('A.mat');
A = single(cell2mat(struct2cell(fid)));

[r,c] = size(A);
delta_S = zeros(62,1); % ��ѳ�
num_F = zeros(62,1);   % �ƻ��������

h = waitbar(0);
for n = 1:62           % J01 ~ J62
    %ը��n����µ��ڽӾ���
    waitbar(n/62,h,['Now processing ' num2str(n) ' / ' num2str(62)]);
    
    A_new = A;
    A_new(n, :) = inf;
    A_new(:, n) = inf;
    A_new(n, n) = 0;
    
    ind = find(A(n,:) ~= inf & A(n,:) ~= 0);
    ind2 = ind(ind <= 70);             % �ҵ���n���ڵĽڵ�ı��
    num_F(n) = length(ind(ind >= 71)); % ���ڵ���71��Ϊ����������
    
    if length(ind2) >= 2
        comb = nchoosek(ind2,2);       % �±��������
        for i = 1:size(comb,1)
            [dist1] = Floyd(A,comb(i,1),comb(i,2));      % ����δը��ǰ��������̾����·��
            [dist2] = Floyd(A_new,comb(i,1),comb(i,2));  % ����ը�ٺ�������µ���̾����·��
            delta_S(n) = delta_S(n) + (dist2 - dist1);   % �����ѳ�
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