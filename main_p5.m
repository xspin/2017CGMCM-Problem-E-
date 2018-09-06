clear
clc

plot_flag = true;

alpha = 1; 
beta = 2;
gamma =100;

pathFile = './J01-62.csv';
posFile = './position.xls';
[pos, ~, raw] = xlsread(posFile);
fid = fopen(pathFile);
t = textscan(fid,'%s', 6*62);
t = t{1};          % Adjacent Table
fclose(fid);
fid = fopen('data/p5.txt', 'w');

% J1-62, D1-2, Z01-06, F01-60
% 1-62, 63-64, 65-70, 71-130
global index_J index_D index_Z index_F;
index_J = [idx('J01'):idx('J62')];
index_D = [idx('D01'):idx('D02')];
index_Z = [idx('Z01'):idx('Z06')];
index_F = [idx('F01'):idx('F60')];

tv = 0;
k = 1000;
while k>0
    fi = randperm(60);
    fi = fi(1:48);
    F_pos = pos(index_F(fi),:);
    vr = variance(F_pos);
    if vr>tv
        tv = vr;
        tfi = fi;
    end
    k = k-1;
end
fprintf(fid, 'variance: %f\n', tv);
index_F = index_F(sort(fi));

global n_J n_D n_Z n_F;
n_J = length(index_J);
n_D = length(index_D);
n_Z = length(index_Z);
n_F = length(index_F);

global v n_ABC;

M = 130;            % # of nodes
A = ones(M)*Inf;    % Adjacent Matrix
for k=1:length(t)
    ss = strsplit(t{k}, ',');
    i = idx(ss{1,1});
    A(i,i) = 0;
    for s = ss(2:end)
        if strcmp(s{1}, '') 
            break
        end
        j = idx(s{1});
        A(i,j) = norm(pos(i,:)-pos(j,:));
        A(j,i) = A(i,j);
    end
end

if plot_flag
    fh = [];
    for k = 1:3
        figure; 
        fh=[fh, axes];
        hold on
%         xlabel('距离(km)');
%         ylabel('距离(km)');
        for i=1:M
            color = 'black';
            if i<63
                f = 'ko'; %J
            elseif i<65
                f = 'b*'; %D
            elseif i<71
                f = 'rd'; %Z
                color = 'r';
            else
                f = 'ms'; %F
            end
            plot(fh(k), pos(i,1), pos(i,2), f)
            text(fh(k), pos(i,1)+1, pos(i,2)+1, raw{i,1}, 'Fontsize', 7, 'Color', color);
            for j=i+1:M
                if A(i,j)<Inf
                    plot(fh(k), [pos(i,1) pos(j,1)], [pos(i,2) pos(j,2)], '-.','Color',[0.5 0.3 0.1])
                end
            end
        end
    end
end

% speeds
v1 = [70, 60, 50];
v = [45, 35, 30];
vr = v ./ v1;

n_ABC = [3,3 ; 3,3; 6,6];

% adjacent matrixes
DA = A;
DB = A;
DC = A;
DABC = {DA, DB, DC};

% J1-J20
for k=1:3
    for i=idx('J1'):idx('J19')
        DABC{k}(i,i+1) = DABC{k}(i,i+1)*vr(k);
        DABC{k}(i+1,i) = DABC{k}(i,i+1);
    end
end
TABC = {};
for k=1:3
    TABC{k} = DABC{k}/v(k);
end

DistABC = {};
TraceABC = {};
DF1 = {};
DF2 = {};
traceZ = {};
tDF1 = {};
tDF2 = {};
I1 = {};
I2 = {};
% loop for computing shortest paths and sort them
for k=1:3 % A B C
    [DistABC{k}, TraceABC{k}] = Floyd(DABC{k});
    [DF1{k}, DF2{k}, traceZ{k}] = spath_DFZF(DistABC{k}, alpha, beta, gamma);
    tDF1{k} = reshape(DF1{k}', [1 n_F*n_F]);
    tDF2{k} = reshape(DF2{k}', [1 n_F*n_F]);
    [tDF1{k}, I1{k}] = sort(tDF1{k});
    [tDF2{k}, I2{k}] = sort(tDF2{k});
end

dist_total = 0;
time_total = 0;
Fe = [];
color = ['r'; 'b'; 'g'];
type = 'ABC';
if plot_flag
    title(fh(1), 'Step 1: D->F');
    title(fh(2), 'Step 2: F->Z');
    title(fh(3), 'Step 3: Z->F');
end

fig = [];
PathABC = {};
f1 = [];
f2 = [];
for k=[3,2,1]
    % greedy algorithm
    [DFZF, ds, Fe] = select(traceZ{k}, tDF1{k}, I1{k}, tDF2{k}, I2{k}, n_ABC(k,:), Fe);
%    DFZFABC{k} = DFZF;
    f1 = [f1; DFZF(:,2)];
    f2 = [f2; DFZF(:,4)];
    for i=1:size(DFZF,1) % travel each path
        p = DFZF(i,:);
        for t=1:length(p)-1
            sp = spath(DistABC{k}, TraceABC{k}, p(t:t+1));
            PathABC{vidx(k, i), t} = sp;
            if plot_flag
                tpos = shift(pos(sp,:), 2, k);
%                 plot(fh(t), tpos(:,1),tpos(:,2), color(k,:), 'LineWidth', 1);
                fig(k,t) = drawArrow(fh(t), tpos(:,1),tpos(:,2), 'Color', color(k,:),'MaxHeadSize', 1, 'LineWidth', 0.5,'AutoScale','off','LineStyle', '-.');
            end
        end
    end
end

F1_pos = pos(f1,:);
F2_pos = pos(f2,:);
F_pos = [F1_pos; F2_pos];
fprintf(fid, 'VarF1:%f, VarF2:%f, VarF:%f\n', variance(F1_pos), variance(F2_pos), variance(F_pos));
for p = F1_pos'
    plot(fh(1), p(1), p(2), 'bp', 'MarkerFaceColor','b');
end
for p = F2_pos'
    plot(fh(3), p(1), p(2), 'rp', 'MarkerFaceColor','r');
end
for p = F_pos'
    plot(fh(2), p(1), p(2), 'kp', 'MarkerFaceColor','k');
end
% varf = @(x) sqrt(sum((x-mean)) / (size(x,1)-1));


% fprintf('\nTotal: d = %.3f km, t = %.3f h\n', dist_total, time_total);

if plot_flag
    lg = legend(fig(:,1), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,2), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,3), 'A', 'B', 'C','Location','southwest');
    saveas(fh(1), 'fig/p5-step1.png');
    saveas(fh(2), 'fig/p5-step2.png');
    saveas(fh(3), 'fig/p5-step3.png');
end

count = zeros(1,n_Z);
for i=1:24
    z = PathABC{i, 3}(1);
    j = find(z==index_Z);
    count(j) = count(j)+1;
end
count

[table, time, extime] = schedule(PathABC, TABC, 0);

fprintf(fid, 'Total Time: %fh\n', time);
fprintf(fid, 'Time:\n');
fprintf(fid, ' F1: %fh\n', max(extime(1,:)));
fprintf(fid, ' F2: %fh\n', max(sum(extime(2:3,:))));
fprintf(fid, ' F1F2: %fh\n', max(sum(extime)));

tb = gettable(table, extime,'data/p5-table.csv', 'data/p5-table.txt');
% time

fclose(fid);



