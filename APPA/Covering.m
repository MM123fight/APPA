% Covering:
% min  c'x 
% s.t. -Ax <= -e, x>= 0;
% A:m-by-n
% Input: m,n,d(density);
rand('seed',sum(100*clock));
m = 4;
n = 3;
d = 1;
tic;
A=sprand(m,n,d);
A=round(A);
A(:, all(A==0,1)) = [];
A(all(A==0,2), :) = [];
[m,n] = size(A);
A = -A;
b = -ones(m,1);
c = rand(n,1);

toc0=toc;
format short e;

mi = m; me = 0;
nb = n; nf = 0;
Ai = A; bi = full(b); 
Ae = sparse(me,n); be = zeros(me,1);
c = full(c);
experiment_name = ['C',num2str(mi),'_', num2str(n), '_' , num2str(d)];
experiment = [pwd, '/data/', experiment_name];
mkdir(experiment);
dataWrite(experiment,nb,nf,mi,me,n,Ae,Ai,be,bi,c);

mkdir([experiment, '/result']);
