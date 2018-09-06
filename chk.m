function [count, wtime, gtime] = chk(TABC, PathABC, table)
fprintf('\n\n****** checking ******\n')
record = cell(130);
wtime=0;
gtime = 0;
for i=1:length(table)
    for t=1:size(table{i},1)-1
        e = table{i}(t:t+1,1);
        record{e(1),e(2)} = [record{e(1),e(2)}; i, table{i}(t,3)*60, table{i}(t+1,2)*60];
        dt = table{i}(t,3) - table{i}(t,2);
        if dt<0
            fprintf('Error:Negative wait time [%.3f] in node %s of index %d\n',dt,idx(e(1)), i);
        end
        if t>1
            wtime = wtime + dt;
        end
        dt = table{i}(t+1,2) - table{i}(t,3);
        if dt<0
            fprintf('Error:Negative run time [%.3f] in edge (%s,%s) of index %d\n',dt,idx(e(1)),idx(e(2)), i);
        end
        gtime = gtime + dt;
    end
    dt = table{i}(end,3) - table{i}(end,2);
    if dt<0
        fprintf('Error:Negative wait time [%.3f] in node %s of index %d(end)\n',dt,idx(e(1)), i);
    end
    wtime = wtime + dt;
%     fprintf('fuck time: %f\n',table{i}(end,3)-table{i}(1,3));
end


count = 0;
delta = 0.0000;
for i=1:size(record,1)
    for j=i+1:size(record,1)
    % check the opsite directions
        for s=1:size(record{i,j},1)
            v1 = record{i,j}(s,1);
            t1 = record{i,j}(s,2:3);
            for t=1:size(record{j,i},1)
                v2 = record{j,i}(t,1);
                t2 = record{j,i}(t,2:3);
                if ((t1(1)<t2(1) && t2(1)<t1(2))||(t1(1)<t2(2) && t2(2)<t1(2))) && v1~=v2
                    count = count +1;
                    fprintf('>< conflict:\n %d & %d in edge (%s, %s)\n', v1, v2, idx(i), idx(j));
                    fprintf('   time:(%.4f,%.4f), (%.4f,%.4f)\n',t1(1), t1(2),t2(1), t2(2));
                end
            end
        end
    % check the same direction of i->j
        for s=1:size(record{i,j},1)
            v1 = record{i,j}(s,1);
            t1 = record{i,j}(s,2:3);
            for t=s+1:size(record{i,j},1)
                v2 = record{i,j}(t,1);
                t2 = record{i,j}(t,2:3);
                if (t2(1)>t1(1) && t2(2)<t1(2))||(t2(1)<t1(1)&&t2(2)>t1(2))
                    count = count +1;
                    fprintf('>> conflict:\n %d & %d in edge (%s, %s)\n', v1, v2, idx(i), idx(j));
                    fprintf('   time:(%.4f,%.4f), (%.4f,%.4f)\n',t1(1), t1(2),t2(1), t2(2));
                end 
            end
        end
    % check the same directions of j->i
        for s=1:size(record{j,i},1)
            v1 = record{j,i}(s,1);
            t1 = record{j,i}(s,2:3);
            for t=s+1:size(record{j,i},1)
                v2 = record{j,i}(t,1);
                t2 = record{j,i}(t,2:3);
                if (t2(1)>t1(1) && t2(2)<t1(2))||(t2(1)<t1(1)&&t2(2)>t1(2))
                    count = count +1;
                    fprintf('>> conflict:\n %d & %d in edge (%s, %s)\n', v1, v2, idx(j), idx(i));
                    fprintf('   time:(%.4f,%.4f), (%.4f,%.4f)\n',t1(1), t1(2),t2(1), t2(2));
                end 
            end
        end
    end
end
