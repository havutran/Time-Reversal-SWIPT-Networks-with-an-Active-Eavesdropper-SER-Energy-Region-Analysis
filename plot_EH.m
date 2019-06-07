clear;
Mt=4;
mu = 0.5;
rhot = 0;
rhol = 0;
PGdb = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18 -22.4 -26.7];
PG   = 1/(40^3)*10.^(PGdb/10);
sum_p = sum(PG);
L=18;
xi = linspace(0,1,11);
pdb = 35;
p = 10.^(pdb/10);
Me = [1 3 5];
pEdb = 7;
pE = 1/(5^3)*10^(pEdb/10).*Me;

[ Pisia ] = ISI_Analytical( p, PGdb, Mt, 0, 0,L );

[ Pisi, Psig ] = ISI_Simulation(  p, PGdb, Mt, 0, 0, L );

[ Pisi_o, Psig_o ] = ISI_Simulation_original (  p, PGdb, Mt, 0, 0, L );

for k = 1: 3
for i = 1: 11
    EHa(k,i) = mu*(1-xi(i))*(Pisia + p*Mt*sum_p + (2*L-1)*pE(k));

    EH(k,i) = mu*(1-xi(i))*(Pisi + p*Mt*sum_p + (2*L-1)*pE(k));

    EHo(k,i) = mu*(1-xi(i))*(Pisi_o + p*Mt*sum_p + (2*L-1)*pE(k));
end
end

%for k = 1: 3
%EHa(k,:) = 10.*log10 (EHa(k,:));
%EH(k,:) = 10.*log10 (EH(k,:));
%EHo(k,:) = 10.*log10 (EHo(k,:));
%end

figure(1); clf;
plot (xi, EHa(1,:),'-k',xi, EH(1,:),'*b',xi, EHo(1,:),'or', xi, EHa(2,:),'-k', xi, EHa(3,:),'-k', xi, EHa(2,:),'*b', xi, EHa(3,:),'*b', xi, EHo(2,:),'or', xi, EHo(3,:),'or','linewidth',1,'MarkerSize',8)
%plot (SNRdb, SINRdb(1,:),'*k',SNRdb, SINRdba(1,:),'-k',SNRdb, SINRdb(2,:),'*k',SNRdb, SINRdba(2,:),'-k',SNRdb, SINRdb(3,:),'*k',SNRdb, SINRdba(3,:),'-k','linewidth',1,'MarkerSize',8)
%grid on;
xlabel('\rho (power splitting factor)')
ylabel('Harvested Energy (mW)')
legend('Analytical - Average harvested energy','Simulation - Average effective harvested energy','Simulation - Average harvested energy',1 )


