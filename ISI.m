%Pisi
%Analytical
clear;
Mt=3;
rhot=0.4;
rhol=0.7;
N = 1;
L = 18;
%PGdb = zeros(1,L);
PGdb =  [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18 -22.4 -26.7];
PG   = 10.^(PGdb/10);
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

Pisia = 2*(w + v + z)/(Mt*sum_p)

%Simulation
P = (PG.^(1/2))'*PG.^(1/2);
x = 1;
psi=0;
S=0;
S1=0;
Sx=0;
Sx1=0;
for k1=1:1000000
%i.i.d channel
hw = manual_channel(L,Mt,1,zeros(1,L));
hw = hw.';

%Correlated channel
h = (RL.*P)^(1/2)*hw*RT^(1/2);

%Equivalent channel
for m=1:Mt
    heq(:,m) = conv(conj(flipud(h(:,m)))/sqrt((Mt*sum_p)), h(:,m));
    h_eq(:,m) = abs(heq(:,m)).^2;
end
hteq = sum(heq,2);
P_heq = abs(hteq).^2;
S=S + P_heq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ht_eq = sum(h_eq,2);
S1 = S1+ ht_eq;

E1 = h(1,1)*conj(h(2,1))/sqrt(PG(1)*(PG(2)));
E = h*h'/(Mt);
Sx = Sx + E;
Sx1 = Sx1 + E1;



%SNR

%SNR = 10.^(SNRdb/10);
%noise1 = P_heq(1,L)./SNR;
%for i=1:length(SNR)
%tempx(1,i) = P_heq(1,L)/(  (sum(P_heq )- P_heq(1,L)  ) + (N-1)*sum(P_heq1) + noise1(1,i));
%end

%Sx = Sx + tempx;
end
S=S/k1;
S1=S1/k1;
Psig = S(L,1);
Pisi = (sum(S) - Psig)

part = sum (S1) - S1(L,1)

Sx=Sx/k1
Sx1=Sx1/k1
















