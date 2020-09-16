% lpgen: Generate random solvable lp:
% min  c'x 
% s.t. Ax =< b; 
% A:m-by-n
% Input: m,n,d(density);

rand('seed',sum(100*clock));
m = 4;
n = 3;
d = 1;
pl=inline('(abs(x)+x)/2');%pl(x) = (|x|+x)/2
tic;
A=sprand(m,n,d);
A=100*(A-0.5*spones(A));
A(:, all(A==0,1)) = [];
A(all(A==0,2), :) = [];
[m,n] = size(A);
% Sparse uniformly distributed random matrix, i.e, every term in [0,1];
% R = spones(A) generates a matrix R with the same sparsity structure as A, 
% but with 1's in the nonzero positions
% A:randomly generated matrix with every term in [-50,50]
u=sparse(10*pl(rand(m,1)-(m-3*n)/m));
x=10*spdiags((sign(pl(rand(n,1)-rand(n,1)))),0,n,n)*(rand(n,1)-rand(n,1));
% dual solution: approximately 3n of which are positvie and m-3n of which
%                are zeros with every term in [0,10]             
% primal solution: approximately half of which are zeros 
%                  with every term in [-10,10];
%c=-A'*u;b=A*x+spdiags((ones(m,1)-sign(pl(u))),0,m,m)*10*ones(m,1);
c=-A'*u;
b=A*x+10*ones(m,1);
%b=A*x+spdiags((ones(m,1)-sign(pl(u))),0,m,m)*10*ones(m,1);
toc0=toc;
format short e;

mi = m; me = 0;
nb = 0; nf = n;
Ai = A; bi = full(b); 
Ae = sparse(me,n); be = zeros(me,1);
c = full(c);
experiment_name = ['Mat',num2str(mi),'_', num2str(n), '_' , num2str(d)];
experiment = [pwd,'/data/',experiment_name];
mkdir(experiment);
dataWrite(experiment,nb,nf,mi,me,n,Ae,Ai,be,bi,c);
mkdir([experiment, '/result']);
