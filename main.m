clear
clc

plot_flag = false;

pathFile = './J01-62.csv';
posFile = './position.xls';
[pos, ~, raw] = xlsread(posFile);
fid = fopen(pathFile);
t = textscan(fid,'%s', 6*62);
t = t{1};          % Adjacent Table
fclose(fid);
% J1-62, D1-2, Z01-06, F01-60
% 1-62, 63-64, 65-70, 71-130
global index_J index_D index_Z index_F;
index_J = [idx('J01'):idx('J62')];
index_D = [idx('D01'):idx('D02')];
index_Z = [idx('Z01'):idx('Z06')];
index_F = [idx('F01'):idx('F60')];
global n_J n_D n_Z n_F;
n_J = length(index_J);
n_D = length(index_D);
n_Z = length(index_Z);
n_F = length(index_F);

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

global v n_ABC;
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
    [DistABC{k}, TraceABC{k}] = Floyd(TABC{k});
    [DF1{k}, DF2{k}, traceZ{k}] = spath_DFZF(DistABC{k}, 1,2,100);
    tDF1{k} = reshape(DF1{k}', [1 n_F*n_F]);
    tDF2{k} = reshape(DF2{k}', [1 n_F*n_F]);
    [tDF1{k}, I1{k}] = sort(tDF1{k});
    [tDF2{k}, I2{k}] = sort(tDF2{k});
end


Time = [0,0,0];

t1 = {[10:12], [7:9], [1:6]};
t2 = {[19:24], [13:18], [1:12]};
for k=1:3
    D12F = [];
    ZF = [];
    for i=index_F
        D12F = [D12F; DistABC{k}(i, index_D)];
        ZF = [ZF, DistABC{k}(i, index_Z)];
    end
    [d1f, I1] = sort(D12F(:,1));
    [d2f, I2] = sort(D12F(:,2));
    [zf, I3] = sort(ZF);
    D1F_12 = d1f(1:12);
    D2F_12 = d2f(1:12);
    zf_24 = zf(1:24);
    Time(k) = Time(k) + sum(D1F_12(t1{k}))+sum(D2F_12(t1{k}))+ sum(zf_24([t1{k}, 12+t1{k}]));
end

Time
sum(Time)
return

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

Z06	356	404.4	J28	461.7	461.7	F19	477.8	456.8

%

for k = [3,2,1]
    % select 12 sortest path from Di -> F, and 24 from F -> Z
end
% selected 




fig = [];
PathABC = {};
for k=[3,2,1]
    % greedy algorithm
    [DFZF, ds, Fe] = seslect(traceZ{k}, tDF1{k}, I1{k}, tDF2{k}, I2{k}, n_ABC(k,:), Fe);
%    DFZFABC{k} = DFZF;
    for i=1:size(DFZF,1) % travel each path
        p = DFZF(i,:);
        for t=1:length(p)-1
            sp = spath(DistABC{k}, TraceABC{k}, p(t:t+1));
            PathABC{vidx(k, i), t} = sp;
            if plot_flag
                tpos = shift(pos(sp,:), 2, k);
%                 plot(fh(t), tpos(:,1),tpos(:,2), color(k,:), 'LineWidth', 1);
%                 annotation(fh(t), 'textarrow',pos(:,1),tpos(:,2),'LineStyle','-','color',[1 0 0]);
                fig(k,t) = drawArrow(fh(t), tpos(:,1),tpos(:,2), 'Color', color(k,:),'MaxHeadSize', 2, 'LineWidth', 1,'AutoScale','off');
            end
        end
    end
    
%     fprintf('Type %s:\n d = %.3f km, t = %.3f h\n', type(k), sum(ds), sum(ds)/v(k))
%     fprintf(' min: d=%.3fkm, t=%.3fh\n', min(ds), min(ds)/v(k));
%     fprintf(' max: d=%.3fkm, t=%.3fh\n', max(ds), max(ds)/v(k));
%     fprintf(' avg: d=%.3fkm, t=%.3fh\n\n', mean(ds), mean(ds)/v(k));
%     dist_total = dist_total + sum(ds);
%     time_total = time_total + sum(ds)/v(k);
    %break;
end

% fprintf('\nTotal: d = %.3f km, t = %.3f h\n', dist_total, time_total);

if plot_flag
    lg = legend(fig(:,1), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,2), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,3), 'A', 'B', 'C','Location','southwest');
    saveas(fh(1), 'fig/p1-step1.png');
    saveas(fh(2), 'fig/p1-step2.png');
    saveas(fh(3), 'fig/p1-step3.png');
    close all
end
% close all

count = zeros(1,n_Z);
for i=1:24
    z = PathABC{i, 3}(1);
    j = find(z==index_Z);
    count(j) = count(j)+1;
end
count
% sum(count)

[table, time, extime] = schedule(PathABC, TABC, 0);
time
sum(sum(extime))
tb = gettable(table, extime,'data/p1-table.csv', 'data/p1-table.txt');
% time





