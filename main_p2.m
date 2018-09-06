clear
clc

plot_flag = true;
loop = true;

pathFile = './J01-62.csv';
posFile = './position.xls';

[pos, ~, raw] = xlsread(posFile);

fid = fopen(pathFile);
tt = textscan(fid,'%s', 6*62);
tt = tt{1};          % Adjacent Table
fclose(fid);

fid=fopen('data/p2.txt','w');
alpha = 1;
beta = 2;
gamma = 100;

% J1-62, D1-2, Z01-06, F01-60
% 1-62, 63-64, 65-70, 71-130

% J25、 J34、 J36、 J42、 J49 
bJ2Z = [idx('J25'), idx('J34'), idx('J36'), idx('J42'), idx('J49')];
% J2Z = J2Z([]);
% J2Z = [];
temp = [];


for jz1=1:length(bJ2Z)
    if ~loop
        break
    end
    for jz2=jz1+1:length(bJ2Z)
        J2Z = bJ2Z([jz1 jz2]);
        fprintf(fid,'%s,%s\n',idx(J2Z(1)),idx(J2Z(2)));
        global index_J index_D index_Z index_F;
        index_J = [idx('J01'):idx('J62')];
        index_D = [idx('D01'):idx('D02')];
        index_Z = [idx('Z01'):idx('Z06')];
        index_F = [idx('F01'):idx('F60')];

        index_Z = [index_Z, J2Z];
        for j = J2Z
            i = find(index_J==j);
            index_J(i) = [];
        end

        global n_J n_D n_Z n_F;
        n_J = length(index_J);
        n_D = length(index_D);;
        n_Z = length(index_Z);;
        n_F = length(index_F);;

        M = 130;            % # of nodes
        A = ones(M)*Inf;    % Adjacent Matrix
        for k=1:length(tt)
            ss = strsplit(tt{k}, ',');
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

        global v;
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
            [DF1{k}, DF2{k}, traceZ{k}] = spath_DFZF(DistABC{k},alpha,beta,gamma);
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

        PathABC = {};
        fig = [];
        for k=[3, 2, 1]
            % greedy algorithm
            [DFZF, ds, Fe] = select(traceZ{k}, tDF1{k}, I1{k}, tDF2{k}, I2{k}, n_ABC(k,:), Fe);
        %    DFZFABC{k} = DFZF;
            for i=1:size(DFZF,1) % travel each path
                p = DFZF(i,:);
                for t=1:length(p)-1
                    sp = spath(DistABC{k}, TraceABC{k}, p(t:t+1));
                    PathABC{vidx(k, i), t} = sp;
                end
            end
        %     fprintf('Type %s:\n d = %.3f km, t = %.3f h\n', type(k), sum(ds), sum(ds)/v(k))
        %     fprintf(' min: d=%.3fkm, t=%.3fh\n', min(ds), min(ds)/v(k));
        %     fprintf(' max: d=%.3fkm, t=%.3fh\n', max(ds), max(ds)/v(k));
        %     fprintf(' avg: d=%.3fkm, t=%.3fh\n\n', mean(ds), mean(ds)/v(k));
            dist_total = dist_total + sum(ds);
            time_total = time_total + sum(ds)/v(k);
            %break;
        end
        fprintf('\nTotal: d = %.3f km, t = %.3f h\n', dist_total, time_total);

        
        count = zeros(1,n_Z);
        number = cell(1,n_Z);
        for i=1:24
            z = PathABC{i, 3}(1);
            j = find(z==index_Z);
            count(j) = count(j)+1;
            if isempty(number{j})
                number{j} = idx(z);
            end
        end
        for i=1:length(count)
            fprintf(fid,' %s : %d\n',number{i}, count(i));
        end
        % % sum(count)
        TABC = {};
        for k=1:3
            TABC{k} = DABC{k}/v(k);
        end
        [~, time] = schedule(PathABC, TABC, 0);
        fprintf(fid, 'Time: %.3f\n',time);
        fprintf(fid, '**********************************\n');
    end
end

%%
% **********************************
% J25,J49
% Total: d = 838.920 km, t = 23.329 h
%  Z01 : 3
%  Z02 : 1
%  Z03 : 1
%  Z04 : 3
%  Z05 : 3
%  Z06 : 5
%  J25 : 5
%  J49 : 3
%  >> dt=0.00
%  >> dt=0.12
% Time: 125.687
if ~plot_flag
    return
end
% J2Z = bJ2Z([jz1 jz2]);
J2Z = [idx('J25'), idx('J49')];
fprintf('%s,%s\n',idx(J2Z(1)),idx(J2Z(2)));
global index_J index_D index_Z index_F;
index_J = [idx('J01'):idx('J62')];
index_D = [idx('D01'):idx('D02')];
index_Z = [idx('Z01'):idx('Z06')];
index_F = [idx('F01'):idx('F60')];

index_Z = [index_Z, J2Z];
for j = J2Z
    i = find(index_J==j);
    index_J(i) = [];
end

global n_J n_D n_Z n_F;
n_J = length(index_J);
n_D = length(index_D);;
n_Z = length(index_Z);;
n_F = length(index_F);;

M = 130;            % # of nodes
A = ones(M)*Inf;    % Adjacent Matrix
for k=1:length(tt)
    ss = strsplit(tt{k}, ',');
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
    [DF1{k}, DF2{k}, traceZ{k}] = spath_DFZF(DistABC{k},alpha,beta,gamma);
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

PathABC = {};
fig = [];
for k=[3, 2, 1]
    % greedy algorithm
    [DFZF, ds, Fe] = select(traceZ{k}, tDF1{k}, I1{k}, tDF2{k}, I2{k}, n_ABC(k,:), Fe);
%    DFZFABC{k} = DFZF;
    for i=1:size(DFZF,1) % travel each path
        p = DFZF(i,:);
        for t=1:length(p)-1
            sp = spath(DistABC{k}, TraceABC{k}, p(t:t+1));
            PathABC{vidx(k, i), t} = sp;
            if plot_flag
                tpos = shift(pos(sp,:), 2, k);
%                         plot(fh(t), tpos(:,1),tpos(:,2), color(k,:), 'LineWidth', 1);
                fig(k,t) = drawArrow(fh(t), tpos(:,1),tpos(:,2), 'Color', color(k,:),'MaxHeadSize', 2, 'LineWidth', 1,'AutoScale','off');
            end
        end
    end
%     fprintf('Type %s:\n d = %.3f km, t = %.3f h\n', type(k), sum(ds), sum(ds)/v(k))
%     fprintf(' min: d=%.3fkm, t=%.3fh\n', min(ds), min(ds)/v(k));
%     fprintf(' max: d=%.3fkm, t=%.3fh\n', max(ds), max(ds)/v(k));
%     fprintf(' avg: d=%.3fkm, t=%.3fh\n\n', mean(ds), mean(ds)/v(k));
    dist_total = dist_total + sum(ds);
    time_total = time_total + sum(ds)/v(k);
    %break;
end
fprintf('\nTotal: d = %.3f km, t = %.3f h\n', dist_total, time_total);

if plot_flag
    lg = legend(fig(:,1), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,2), 'A', 'B', 'C','Location','southwest');
    legend(fig(:,3), 'A', 'B', 'C','Location','southwest');
    saveas(fh(1), 'fig/p2-step1.png');
    saveas(fh(2), 'fig/p2-step2.png');
    saveas(fh(3), 'fig/p2-step3.png');
    close all
end

count = zeros(1,n_Z);
number = cell(1,n_Z);
for i=1:24
    z = PathABC{i, 3}(1);
    j = find(z==index_Z);
    count(j) = count(j)+1;
    if isempty(number{j})
        number{j} = idx(z);
    end
end
for i=1:length(count)
    fprintf(fid,' %s : %d\n',number{i}, count(i));
end
% % sum(count)
TABC = {};
for k=1:3
    TABC{k} = DABC{k}/v(k);
end
[~, time] = schedule(PathABC, TABC, 0);
fprintf(fid,'Time: %.3fh\n',time);
fprintf(fid,'**********************************\n');

[table, time, extime] = schedule(PathABC, TABC, 0);
tb = gettable(table, extime, 'data/p2-table.csv', 'data/p2-table.txt');
fclose(fid);

close all

