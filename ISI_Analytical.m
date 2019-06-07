function [ Pisia ] = ISI_Analytical( p, PGdb, Mt, rhot, rhol,L )

PG   = (1/(10^3))*10.^(PGdb/10);
sum_p = sum(PG);
RT = CorrMatrix_interclass(Mt,rhot);
RL = CorrMatrix_interclass(L,rhol);
%RL = RL.*P;
w=0;
v=0;
z=0;

for k = 1: L-1
    for m1=1:Mt
        for m2=1:Mt
            if m1~=m2
                for l1=1:k
                    for l2=1:k 
                        w = w + RL(k+1-l1,L+1-l1)*RL(k+1-l2,L+1-l2)*sqrt(PG(k+1-l1)*PG(L+1-l1))*sqrt(PG(k+1-l2)*PG(L+1-l2)) + RT(m1,m2)*RL(k+1-l1,k+1-l2)*RT(m1,m2)*RL(L+1-l1,L+1-l2)*sqrt(PG(k+1-l1)*PG(L+1-l1))*sqrt(PG(k+1-l2)*PG(L+1-l2));
                    end
                end
            end
        end
    end
    for l3=1:k
       v = v + Mt*(RL(k+1-l3,L+1-l3)^2+1)*PG(k+1-l3)*PG(L+1-l3); 
    end
    for l4=1:k
        for l5=1:k
            if l4 ~= l5
                z = z + Mt*RL(k+1-l4,L+1-l4)*RL(k+1-l5,L+1-l5)*sqrt(PG(k+1-l4)*PG(L+1-l4))*sqrt(PG(k+1-l5)*PG(L+1-l5)) + Mt*RL(k+1-l4,k+1-l5)*RL(L+1-l4,L+1-l5)*sqrt(PG(k+1-l4)*PG(L+1-l4))*sqrt(PG(k+1-l5)*PG(L+1-l5));
            end
        end
    end
end

Pisia = 2*p*(w + v + z)/(Mt*sum_p);

end

