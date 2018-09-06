clear
clc

plot_flag = true;

pathFile = './J01-62.csv';
posFile = './position.xls';
[pos, ~, raw] = xlsread(posFile);
fid = fopen(pathFile);
t = textscan(fid,'%s', 6*62);
t = t{1};          % Adjacent Table
fclose(fid);

fid = fopen('data/p3.txt','w');

alpha = 1;
beta = 2;
gamma = 100;
% J1-62, D1-2, Z01-06, F01-60
% 1-62, 63-64, 65-70, 71-130

global index_J index_D index_Z index_F;
index_J = [idx('J01'):idx('J62')];
index_D = [idx('D01'):idx('D02')];
index_Z = [idx('Z01'):idx('Z06')];
index_F = [idx('F01'):idx('F60')];
global n_J n_D n_Z n_F;
n_J = length(index_J);
n_D = length(index_D);;
n_Z = length(index_Z);;
n_F = length(index_F);;

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
                f = 'gs'; %F
            end
            plot(fh(k), pos(i,1), pos(i,2), f)
            text(fh(k), pos(i,1)+1, pos(i,2)+1, raw{i,1}, 'Fontsize', 7, 'Color', color);
            for j=i+1:M
                if A(i,j)<Inf
                    plot(fh(k), [pos(i,1) pos(j,1)], [pos(i,2) pos(j,2)], 'k-.')
                end
            end
        end
    end
end

global v;
% speeds
v1 = [70, 60, 50];
v = [45, 35, 30];
vr = v ./ v1;

% n_ABC = [3,3 ; 3,3; 6,6];

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

color = ['r'; 'b'; 'g'];
type = 'ABC';
if plot_flag
    title(fh(1), 'Step 1: D->F');
    title(fh(2), 'Step 2: F->Z');
    title(fh(3), 'Step 3: Z->F');
end
global n_ABC;
n_ABC = [3,3 ; 3,3; 6,6];
n_C = [0,3  ; 1,2 ; 3,0; 2,1];

temp = [];
for i_nc=1:size(n_C, 1)
    tn_ABC = n_ABC;
    tn_ABC(3,:) = n_ABC(3,:) - n_C(i_nc,:);
    PathABC = {};
    Fe = [];
    for k=[3, 2, 1]
        % greedy algorithm
        [DFZF, ds, Fe] = select(traceZ{k}, tDF1{k}, I1{k}, tDF2{k}, I2{k}, tn_ABC(k,:), Fe);
    %    DFZFABC{k} = DFZF;
        
    
        for i=1:size(DFZF,1) % travel each path
            p = DFZF(i,:);
            for t=1:length(p)-1
                sp = spath(DistABC{k}, TraceABC{k}, p(t:t+1));
                PathABC{vidx(k, i), t} = sp;
                if plot_flag
                    tpos = shift(pos(sp,:), 2, k);                end
            end
        end
%         fprintf('===%d, %d\n',k, size(DFZF,1))
%         fprintf('Type %s:\n d = %.3f km, t = %.3f h\n', type(k), sum(ds), sum(ds)/v(k))
%         fprintf(' min: d=%.3fkm, t=%.3fh\n', min(ds), min(ds)/v(k));
%         fprintf(' max: d=%.3fkm, t=%.3fh\n', max(ds), max(ds)/v(k));
%         fprintf(' avg: d=%.3fkm, t=%.3fh\n\n', mean(ds), mean(ds)/v(k));
%         dist_total = dist_total + sum(ds);
%         time_total = time_total + sum(ds)/v(k);
        %break;
    end
%     fprintf('\nTotal: d = %.3f km, t = %.3f h\n', dist_total, time_total);

    % n_ABC(3,:) from D -> F
    [CDF1, Fe] = spath_XF(DistABC{3}, index_D(1), n_C(i_nc,1), Fe);
    [CDF2, Fe] = spath_XF(DistABC{3}, index_D(2), n_C(i_nc,2), Fe);
    ki = 22;
%     size(PathABC)
    for df = [CDF1', CDF2']
        sp = spath(DistABC{3}, TraceABC{3}, df);
        PathABC{ki, 1} = sp;
        ki = ki+1;
    end
    % J04、 J06、 J08、 J13、 J14、 J15 to F
    ks = [idx('J04'), idx('J06'), idx('J08'), idx('J13'), idx('J14'), idx('J15')];
    [CJF, Fe] = spath_XF(DistABC{3}, ks, 3, Fe);
    ki = 22;
    for jf=CJF'
        sp = spath(DistABC{3}, TraceABC{3}, df);
        PathABC{ki, 2} = [sp(1)];
        PathABC{ki, 3} = sp;
        ki = ki+1;
    end
 
    
    [~, time] = schedule(PathABC, TABC, 0);
    temp = [temp; n_C(i_nc,:), CJF(:,1)', time];
    tp = temp(end,:);
    fprintf(fid, ' D1:%d, D2:%d\n', tp(1), tp(2))
    fprintf(fid, ' %s, %s, %s\n',idx(tp(3)), idx(tp(4)), idx(tp(5)));
    fprintf(fid, ' Time: %.3fh\n',tp(end))
%     [~, time] = schedule(PathABC, TABC, 0)
    fprintf(fid, '---------------\n')
end

[t, i] = min(temp(:,end));
tp = temp(i,:);
fprintf(fid,' Optimal:\n')
fprintf(fid,'  D1:%d, D2:%d\n', tp(1), tp(2))

fprintf(fid,'  %s, %s, %s\n',idx(tp(3)), idx(tp(4)), idx(tp(5)));
fprintf(fid,'  Time: %.3fh\n',tp(end))

fclose(fid);
%J04、 J06、 J08、 J13、 J14、 J15

% table = schedule(PathABC, TABC, 0);
%%


tn_ABC = n_ABC;
tn_ABC(3,:) = n_ABC(3,:) - tp(1:2);
PathABC = {};
Fe = [];
fig=[];
for k=[3, 2, 1]
    % greedy algorithm
    [DFZF, ds, Fe] = select(traceZ{k}, tDF1{k}, I1{k}, tDF2{k}, I2{k}, tn_ABC(k,:), Fe);
    for i=1:size(DFZF,1) % travel each path
        p = DFZF(i,:);
        for t=1:length(p)-1
            sp = spath(DistABC{k}, TraceABC{k}, p(t:t+1));
            PathABC{vidx(k, i), t} = sp;
            if plot_flag
                tpos = shift(pos(sp,:), 2, k);
                fig(k,t) = drawArrow(fh(t), tpos(:,1),tpos(:,2), 'Color', color(k,:),'MaxHeadSize', 2, 'LineWidth', 1,'AutoScale','off');
%                 plot(fh(t), tpos(:,1),tpos(:,2), color(k,:), 'LineWidth', 1);
            end
        end
    end
%         fprintf('===%d, %d\n',k, size(DFZF,1))
%         fprintf('Type %s:\n d = %.3f km, t = %.3f h\n', type(k), sum(ds), sum(ds)/v(k))
%         fprintf(' min: d=%.3fkm, t=%.3fh\n', min(ds), min(ds)/v(k));
%         fprintf(' max: d=%.3fkm, t=%.3fh\n', max(ds), max(ds)/v(k));
%         fprintf(' avg: d=%.3fkm, t=%.3fh\n\n', mean(ds), mean(ds)/v(k));
%         dist_total = dist_total + sum(ds);
%         time_total = time_total + sum(ds)/v(k);
    %break;
end
%     fprintf('\nTotal: d = %.3f km, t = %.3f h\n', dist_total, time_total);
if plot_flag
    lg = legend(fig(:,1), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,2), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,3), 'A', 'B', 'C','Location','southwest');
    saveas(fh(1), 'fig/p3-step1.png');
    saveas(fh(2), 'fig/p3-step2.png');
    saveas(fh(3), 'fig/p3-step3.png');
%     close all
end

% n_ABC(3,:) from D -> F
[CDF1, Fe] = spath_XF(DistABC{3}, index_D(1), n_C(i_nc,1), Fe);
[CDF2, Fe] = spath_XF(DistABC{3}, index_D(2), n_C(i_nc,2), Fe);
ki = 22;
%     size(PathABC)
for df = [CDF1', CDF2']
    sp = spath(DistABC{3}, TraceABC{3}, df);
    PathABC{ki, 1} = sp;
    ki = ki+1;
end
% J04、 J06、 J08、 J13、 J14、 J15 to F
ks = [idx('J04'), idx('J06'), idx('J08'), idx('J13'), idx('J14'), idx('J15')];
[CJF, Fe] = spath_XF(DistABC{3}, ks, 3, Fe);
ki = 22;
for jf=CJF'
    sp = spath(DistABC{3}, TraceABC{3}, df);
    PathABC{ki, 2} = [sp(1)];
    PathABC{ki, 3} = sp;
    ki = ki+1;
end



[table, time, extime] = schedule(PathABC, TABC, 0);
time

tb = gettable(table, extime, 'data/p3-table.csv', 'data/p3-table.txt');





