function [ Pisi, Psig ] = ISI_Simulation_original (  p, PGdb, Mt, rhot, rhol, L )

PG   = (1/(10^3))*10.^(PGdb/10);
sum_p = sum(PG);
RT = CorrMatrix_interclass(Mt,rhot);
RL = CorrMatrix_interclass(L,rhol);
P = (PG.^(1/2))'*PG.^(1/2);
x = 1;
psi=0;
S=0;
S1=0;
Sx=0;
Sx1=0;
for k1=1:100000
%i.i.d channel
hw = manual_channel(L,Mt,1,zeros(1,L));
hw = hw.';

%Correlated channel
h = (RL.*P)^(1/2)*hw*RT^(1/2);

%Equivalent channel
for m=1:Mt
    heq(:,m) = conv(conj(flipud(h(:,m)))/norm(h,'fro'), h(:,m));
    h_eq(:,m) = abs(heq(:,m)).^2;
end
hteq = sum(heq,2);
P_heq = abs(hteq).^2;
S=S + P_heq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ht_eq = sum(h_eq,2);
S1 = S1+ ht_eq;


end
S=S/k1;
S1=S1/k1;
Psig = S(L,1);
Pisi = p*(sum(S) - Psig);

end

