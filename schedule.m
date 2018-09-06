function [table, alltime, exptime] = schedule(PathABC, TABC, flag)
% PathABC{k, i, t}, k=A/B/C, i=1:6/12, t=DF/FZ/ZF
% table 发射车 S1 S2 S3 时间  待机地域编号	出发时刻	道路节点编号	到达时刻	离开时刻	道路节点编号
table = {};
% {k, i}
% A1-A6, B1-B6, C1-C12
% index = [1, 4, 7, 10, 13, 19]
type = 'ABC';
% n = [6  6 12];
v = [45, 35, 30];
ExpTime = [];
global n_Z index_Z;

for k=1:sum(size(PathABC,1))
    table{k} = [];
end
       
StepTime = [];
EndTime = [0, 0, 0];
exptime = zeros(3, size(PathABC,1));
waittime = 0;

% dists = {[],[],[]};
% dist = 0;


%% Step 1: D->F
step = 1;
if flag
    fprintf('\n===================================\n');
    fprintf('\nStep %d: \n',step);
end

dists = [];
for i=1:size(PathABC,1)
    dist = pathdist(TABC, PathABC{i,step}, vidx(i));
    dists = [dists, dist];
end
[~, I] = sort(dists);
ExpTime = 0;
maxtime = max(dists);
for i = I(end:-1:1)
    stpath = PathABC{i, step};
    time = maxtime-dists(i);
    table{i} = [table{i}; double(stpath(1)), 0, time];
    for t=1:length(stpath)-1
        edge = stpath(t:t+1);
        time = time + TABC{vidx(i)}(edge(1), edge(2));
        table{i} = [table{i}; double(edge(2)), time, time];
        [wt, ft] = wait(table, i);
        if ft>0
            fprintf(' speed up: %f\n',ft);
        end
        if wt>0
            table{i}(end-1, 3) = table{i}(end-1, 3)+wt;
            table{i}(end, 2) = table{i}(end, 2)+wt;
            table{i}(end, 3) = table{i}(end, 3)+wt;
            waittime = waittime + wt;
            fprintf(' wait: %f\n', wt);
        elseif wt<0
            fprintf('error: %f',wt);
        end
    end
    if table{i}(end,2)>EndTime(step)
        EndTime(step) = table{i}(end,2);
    end
    ExpTime = ExpTime + table{i}(end, 2)-table{i}(1, 3);
    exptime(step, i) = exptime(step, i)+table{i}(end, 2)-table{i}(1, 3);
end


if flag
    fprintf('Step Time: %.3fh, End Time: %fh, max: %f\n', ExpTime, EndTime(step), maxtime)
end
StepTime = [StepTime, ExpTime];
% table{1:4}


%% Step 2: F->Z
step = 2;
if flag
    fprintf('\n===================================\n');
    fprintf('\nStep 2: \n', step);
end

dists = [];
for i=1:size(PathABC,1)
    dist = pathdist(TABC, PathABC{i,step}, vidx(i))+pathdist(TABC, PathABC{i,step+1}, vidx(i));
    dists = [dists, dist];
end
[~, I1] = sort(dists);
idx_m = I1(end);
endtime = EndTime(step-1) + max(dists)+1/6;
% endtime = endtime
dists = [];
for i=1:size(PathABC,1)
    dist = pathdist(TABC, PathABC{i,step}, vidx(i));
    dists = [dists, dist];
end
[~, I] = sort(dists);
% I(find(I==idx_m)) = [];

queue = cell(1,n_Z);
timeinz = cell(1,n_Z);
timeoutz = cell(1,n_Z);
zi = zeros(length(I));
wtime = 0;
ExpTime = 0;
lefttime = zeros(1,length(I));
delta = 0.00001;
savetime = 0;
length(I1)
length(I)
% for i = [idx_m, I]
for i = I(end:-1:1)
% for i = I
    stpath = PathABC{i, step};
    time = EndTime(step-1);
    table{i}(end, 3) = EndTime(step-1);
    zi(i) = size(table{i},1);
%     table{i} = [table{i}; double(stpath(1)), 0, time];
    for t=1:length(stpath)-1
        edge = stpath(t:t+1);

        time = time + TABC{vidx(i)}(edge(1), edge(2));
        table{i} = [table{i}; double(edge(2)), time, time];
        wt = wait(table, i);
        if wt>0
            table{i}(end-1, 3) = table{i}(end-1, 3)+wt;
            table{i}(end, 2) = table{i}(end, 2)+wt;
            table{i}(end, 3) = table{i}(end, 3)+wt;
            waittime = waittime + wt;
        elseif wt<0
            fprintf('FFFFFFFF: %f',wt);
        end
    end
    if length(stpath)<=1
        tm = endtime - pathdist(TABC, PathABC{i,3}, vidx(i));
        table{i} = [table{i}; PathABC{i, step}(1), tm, tm];
        continue;
    end
    if table{i}(end,2)>EndTime(step)
        EndTime(step) = table{i}(end,2);
    end
     ExpTime = ExpTime + table{i}(end, 2)-EndTime(step-1);
     exptime(step, i) = exptime(step, i) + table{i}(end, 2)-EndTime(step-1);
     
    % loading and waiting or hiding in Z
    table{i}(end,3) = table{i}(end,3)+1/6;
    lefttime(i) = endtime-table{i}(end,3)-pathdist(TABC, PathABC{i,3}, vidx(i));
    z = PathABC{i, 3}(1);
    j = find(z==index_Z);
    if(lefttime(i)<0)
        fprintf('no left time of %d : %f\n',i,lefttime(i));
    end
    timeinz{j} = [timeinz{j}, table{i}(end,2)];

    if length(queue{j})==2
        k1 = queue{j}(1);
        k2 = queue{j}(2);
%         if table{k2}(end,3)<=table{i}(end,2)        %i在k2装载后到
%             if lefttime(k1)>=table{i}(end,2)-table{k1}(end,3)+delta
%                 table{k1}(end,3) = table{i}(end,2);            % 如果k1剩余时间多就等i
%             else
%                 table{k1}(end,3) = table{k1}(end,3)+lefttime(k1); %否则最多到剩余时间处走
%             end
%         else                                    %i在k2装载中到
%             if lefttime(k1)>=table{k1}(end,3)-table{i}(end,3)+delta %如果k1的剩余时间够等i
%                 table{k1}(end,3) = table{i}(end,2);
%             else                                            %否则
%                 table{k1}(end,3) = table{k1}(end,3)+lefttime(k1);
%             end
%         end
        dt = table{i}(end,2)-table{k1}(end,3); % i到时与k1最快走时的时间差
        if lefttime(i)>lefttime(k1)-max(0, -dt)
            % 如果i的剩余时间更多，则i留下，k1走
            table{k1}(end,3) = table{k1}(end,3) + max(0, dt);
            lefttime(k1) = lefttime(k1) - max(0, dt);
            table{i}(end,3) = table{i}(end,3) + max(0, -dt);
            lefttime(i) = lefttime(i) - max(0, -dt);
            exptime(step, i) = exptime(step, i) + max(0, -dt);
            ExpTime = ExpTime + max(0,-dt);
            waittime = waittime + max(0,-dt);
            savetime = savetime - max(0,-dt);
            fprintf('%d wait in Z%d: %f\n',i, j, max(0,-dt));
            if lefttime(i)>lefttime(k2)-max(0, table{k2}(end,3)-table{i}(end,2))
                queue{j}= [k2, i];
            else
                queue{j}= [i, k2];
            end
        else
            % i剩余时间更少，则i走，k1等10min
            exptime(step, k1) = exptime(step, k1) + 1/6;
            ExpTime = ExpTime + 1/6;
            savetime = savetime - 1/6;
%             table{k1}(end,3) = table{i}(end,3)
        end
%         fprintf('%d 2: %d->%d->%d\n',j, i,k2,k1);
    elseif length(queue{j})==1
        k1 = queue{j}(1);
        if table{k1}(end,3)<=table{i}(end,2) % i在k1装载后
            if lefttime(k1)+delta<=table{i}(end,2)-table{k1}(end,3)
                table{k1}(end,3) = table{k1}(end,3)+lefttime(k1);
                queue{j}=[i];
            else
                if lefttime(k1)-(table{i}(end,2)-table{k1}(end,3))>lefttime(i)
                    queue{j} = [i, k1];
                else
                    queue{j} = [k1, i];
                end
            end
        else
            table{i}(end,3) = table{i}(end,3) + table{k1}(end,3)-table{i}(end,2); %i要等k1完成
            if lefttime(k1)>lefttime(i)-(table{k1}(end,3)-table{i}(end,2));
                queue{j} = [i, k1];
            else
                queue{j} = [k1, i];
            end
        end
%         fprintf('%d 1: %d->%d\n',j, i,k1);
    else
        if lefttime(i)-delta>0
            queue{j}=[i];
%             fprintf('%d 0: %d, %.2f->\n',j, i, lefttime(i));
        end
    end
    timeoutz{j} = [timeoutz{j}, table{i}(end,3)];
end

figure
hold on
xlabel('t');
for j=1:n_Z
    fprintf('Z%d: ', j)
    plot(timeinz{j}*60, j*ones(1,length(timeinz{j})), '-*');
    plot(timeoutz{j}*60, j*ones(1,length(timeoutz{j}))+0.1, '-o');
    disp(timeinz{j})
    for i = queue{j}
        lefttime = endtime-table{i}(end,3)-pathdist(TABC, PathABC{i,3}, vidx(i));
        
        if lefttime>0
            table{i}(end,3) = table{i}(end,3)+lefttime;
        elseif lefttime<=0
            fprintf('error: more left time: %f\n',lefttime);
        end
        timeoutz{j}(end) = table{i}(end,3);
%         lefttime = endtime-table{i}(end,3)-pathdist(TABC, PathABC{i,3}, vidx(i))
    end
end
% return
for i=1:length(table)
    savetime = savetime + table{i}(end,3)-table{i}(end,2);
%     lefttime = endtime-table{i}(end,3)-pathdist(TABC, PathABC{i,3}, vidx(i));
%     fprintf('%2d left time: %f\n', i, lefttime)
end

StepTime = [StepTime, ExpTime];

if flag
    fprintf('Step Time: %.3fh, End Time: %fh, Save Time: %fh\n', ExpTime, EndTime(step),savetime)
end
% endtime

%% Step 3: Z->F
step = 3;
if flag
    fprintf('\n===================================\n');
    fprintf('\nStep %d: \n', step);
end
dists = [];
for i=1:size(PathABC,1)
    dist = pathdist(TABC, PathABC{i,step}, vidx(i));
    dists = [dists, dist];
end
[~, I] = sort(dists);
ExpTime = 0;
temp = [];
wtime = [];
endtime2=0;
for i = I(end:-1:1)
    stpath = PathABC{i, step};
    time = table{i}(end, 3);
    temp = time;
%     if table{i}(end,3)<table{i}(end,2)+1/6-0.1
%         fprintf('ERRORRRRRR!')
%         table{i}(end,:)
%     end
%     table{i} = [table{i}; double(stpath(1)), 0, time];
    for t=1:length(stpath)-1
        edge = stpath(t:t+1);
        time = time + TABC{vidx(i)}(edge(1), edge(2));
        table{i} = [table{i}; double(edge(2)), time, time];
        wt = wait(table, i);
        if wt>0
            table{i}(end-1, 3) = table{i}(end-1, 3)+wt;
            table{i}(end, 2) = table{i}(end, 2)+wt;
            table{i}(end, 3) = table{i}(end, 3)+wt;
            waittime = waittime + wt;
        elseif wt<0
            fprintf('FFFFFFFF: %f',wt);
        end
    end
    if table{i}(end,2)>EndTime(step)
        EndTime(step) = table{i}(end,2);
    end
    if table{i}(end, 2)>endtime2
        endtime2 = table{i}(end, 2);
    end
    ExpTime = ExpTime + table{i}(end, 2)-temp;
    exptime(step,i) = exptime(step, i) + endtime-temp;
    table{i}(end,:);
end
% endtime
for i = I(end:-1:1)
    wtime = [wtime, endtime2-table{i}(end,2)];
    table{i}(end, 3) = endtime2;
end
ExpTime = ExpTime+sum(wtime);

StepTime = [StepTime, ExpTime];

if flag
    fprintf('Step Time: %.3fh, End Time: %fh, max: %f\n', ExpTime, EndTime(step), endtime)
end

%% 
if flag
    fprintf('\n===================================\n');
    fprintf('\nTotal Time: %.3fh\n', sum(StepTime))

end

alltime = sum(StepTime);

fprintf('endwait time:\n')
% disp(wtime)
fprintf(' total: %.3f\nallwait time: %.3f\n', sum(wtime), waittime)
fprintf('end time: %.3f, real end time: %.3f\n',endtime, endtime2);
fprintf('step time: %.3f. %.3f, %.3f, total: %.3f\n',StepTime, alltime);

[count, wtime, gtime] = chk(TABC, PathABC, table);
fprintf(' conflict: %d\n wtime: %f, gtime: %f, total:%f\n Exp:%f\n', count, wtime, gtime, wtime+gtime,wtime+gtime-savetime);
    