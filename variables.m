%edge (i, j)
pos(i,:)    % position of node i
A(i, j)     % Adjacent Matrix with distance
type='ABC'
v(k)        % velocity of type k car
DABC{k}     % Adjacent Matrix of type k
DistABC{k}  % shortest path distance
TraceABC{k} % traces of the shortest path
DF1{k}      % shortest path of DFZF from D1
DF2{k}
traceZ{k}   % corresponding Z in DFi{k}
tDFi{k}     % sorted DFi{k}
Ii{k}       % corresponding index in original unsorted tDFi
PathABC{ki,s} % selected paths of type ki in step s
table{ki}  % vehicle ki: [T1 T2 T3; Pi ti1 ti2; ...]
record{i, j} % traces in edge (i->j): vehk, timek; ...
ecf{i,j}    % f r1 r2; ...




[DFZF, ds, Fe] = select(traceZ, tDF1, I1, tDF2, I2, n, Fe)
% DFZF: paths 
% ds: distances
% Fe: selected F nodes
[table] = schedule(PathABC, DABC)
% PathABC{k, i, t}, k=A/B/C, i=1:6/12, t=DF/FZ/ZF
% table{k,i}: [t0, t11, t12, t21, t22, ...,]
[ecf, count]= conchk(record, DABC)
% check if any conflicted
% ecf{i,j} = [f r1 r2 ; ... ]
