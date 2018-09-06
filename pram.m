function  PathABC = pram(DABC, alpha, beta, gamma)

global n_J n_D n_Z n_F n_ABC;
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
    [DF1{k}, DF2{k}, traceZ{k}] = spath_DFZF(DistABC{k},alpha, beta, gamma);
    tDF1{k} = reshape(DF1{k}', [1 n_F*n_F]);
    tDF2{k} = reshape(DF2{k}', [1 n_F*n_F]);
    [tDF1{k}, I1{k}] = sort(tDF1{k});
    [tDF2{k}, I2{k}] = sort(tDF2{k});
end

dist_total = 0;
time_total = 0;
Fe = [];


PathABC = {};
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
end