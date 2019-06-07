function [SER] = SEP_TRP2P_analytical (rhot,rhol,PGdb,xi, Me, pEdb, pAdb)

M=8;
n_r = 1;
n_t = 4;
L=18;
I = eye(18*n_t);
d = 40;
PG   = 10.^(PGdb/10);
theta = linspace (0,pi - pi/M,10000);

sum_p = sum(PG);

pE = Me*10.^(pEdb/10);

pA = 10.^(pAdb/10);

No = (1/(d^3))*pA/(10^(30/10));

%No = 10^(-11);

R_t = CorrMatrix_interclass(n_t,rhot);
R_L = CorrMatrix_interclass(L,rhol);


T = kron(sqrt(PG')*sqrt(PG).*R_L,R_t);

gamma = (1/(d^3))*xi*pA./((1/(4^3))*xi*pE + No);

for k=1:length(gamma)
fun = 0;   
for i=2:length(theta)
    fun = fun +(1/pi)*det(I + gamma(1,k)*(sin(pi/M))^2/(sin(theta(1,i)))^2*T)^(-1);
end

SER(1,k) = fun*(pi - pi/M)/10000;

end
end

