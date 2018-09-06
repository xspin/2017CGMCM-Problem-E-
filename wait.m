function [wt, ft]= wait(table, i)
wt = 0;
ft = 0;
for j=1:length(table)
    if i==j
        continue;
    end
    for k=1:size(table{j},1)-1
        dt = [];
        ti1 = table{i}(end-1, 3);
        ti2 = table{i}(end, 2);
        if table{i}(end-1,1) == table{j}(k,1) && table{i}(end,1) == table{j}(k+1,1)
            tj1 = table{j}(k,3);
            tj2 = table{j}(k+1,2);
            if tj1>ti1 && tj2<ti2
                dt = (tj1-ti1);
            elseif tj1<ti1 && tj2>ti2
                dt = tj2-ti2;
                ft = ft + ti1-tj1;
            end
            pre = sprintf('%s>>%s:', idx(table{j}(k,1)), idx(table{j}(k+1,1)));
%             if vidx(i)<vidx(j)
%                 ft = ft + tj1
%             end
        elseif table{i}(end-1,1) == table{j}(k+1,1) && table{i}(end,1) == table{j}(k,1)
            tj1 = table{j}(k,3);
            tj2 = table{j}(k+1,2);
%             if (tj1>ti1&&tj1<ti2) || (ti1>tj1&&ti1<tj2)
            if (ti1<tj1&&tj1<ti2) || (ti1<tj2&&tj2<ti2) 
                dt = tj2-ti1;
            end
            pre = sprintf('%s><%s:', idx(table{j}(k,1)), idx(table{j}(k+1,1)));
        end
        
        if dt>0
            fprintf(' %s dt=%f\n',pre, dt);
            wt = wt + dt;
            table{i}(end-1, 3) = table{i}(end-1, 3)+dt;
            table{i}(end, 2) = table{i}(end, 2)+dt;
        elseif dt<=0
            fprintf('Warnning %s dt=%f\n',pre, dt);
        end
    end
end