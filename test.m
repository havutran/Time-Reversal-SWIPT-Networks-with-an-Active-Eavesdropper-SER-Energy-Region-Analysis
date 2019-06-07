clear;

Mt=4;
rho1=0;
rho2=0;
N = 1;
L = 6;
PGdb = [0 -1 -9 -10 -15 -20];

PG = 10.^(PGdb/10);
sum_p = sum(PG);

S=0;
S1=0;
Sx=0;

%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_l = CorrMatrix_interclass(N,rho2);

%Simlation
F = 0;
for k=1:100000
%i.i.d channel
hw = manual_channel(6,Mt,1,PGdb);
h = R_t^(1/2)*hw*R_l^(1/2);

F = F + norm(h,'fro')^2;

end
F =F/100000

T = kron(diag(PG),R_t);
trace(T)

trace(diag(PG))*Mt