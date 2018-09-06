function tb = gettable(table, extime,fname, f2)

global n_ABC;
fprintf('Save to %s, %s\n', fname, f2);

fid = fopen(fname,'w');
fid2 = fopen(f2,'w');
n = sum(n_ABC');
type='ABC';
head ={'发射车编号';'待机地域编号';'出发时刻';	'道路节点编号';'到达时刻';'离开时刻'};
fprintf(fid, '%s,%s,%s',head{1},head{2},head{3});
for i=1:7
    fprintf(fid, ',%s,%s,%s',head{4},head{5},head{6});
end
fprintf(fid, '......\n');
delta=60;
extime = extime*delta;
i_veh = 1;
tb = {};
for k=1:3
    di = [];
    for j=1:sum(n(k))
        i = sum(n(1:k-1))+j;
        di = [di, table{i}(1,1)];
    end
    [~, I] = sort(di);
    for j=1:sum(n(k))
        i = sum(n(1:k-1))+I(j);
        %vehcile N.O.
        fprintf(fid, '%s%02d', type(k), j);
        fprintf(fid, ',%s', idx(table{i}(1,1)));
        tb{i_veh} = sprintf('%s%02d', type(k), j);
        tb{i_veh} = [tb{i_veh}, sprintf(':%s', idx(table{i}(1,1)))];
        fprintf(fid, ',%.1f', table{i}(1,3)*delta);
        for r=table{i}(2:end,:)'
            fprintf(fid, ',%s,%.1f,%.1f',idx(r(1)),r(2)*delta,r(3)*delta);
            tb{i_veh} = [tb{i_veh}, sprintf('->%s', idx(r(1)))];
        end
        fprintf(fid, '\n');
        fprintf(fid2, '%s\n', tb{i_veh});
%         fprintf(fid2, '%s|%.1f|%.1f|%.1f||%.1f\n', tb{i_veh},extime(1,i_veh),extime(2,i_veh),extime(3,i_veh),sum(extime(:,i_veh)));
        i_veh = i_veh+1;
    end
end
fclose(fid);
fclose(fid2);

% 1st
% 1,2,100  134.8057
% 
%
% 2nd
%
%
% third
% 