function [SER] = SEP_TRP2P (rhot,rhol,PGdb,xi, Me, pEdb, pAdb)
M = 8;
K = log2(M);
Nt = 1;
L = 18;
PG   = 10.^(PGdb/10);
P = (PG.^(1/2))'*PG.^(1/2);

d = 40;

No = (1/(d^3))*10.^(pAdb/10)/(10^(30/10));

pE = Me*(1/(4^3))*xi*10.^(pEdb/10);

pA = (1/(d^3))*xi*10.^(pAdb/10);

switch M
case 2                    % BPSK symbols
   sym_map=[1;-1];
   Es = pA;
case 4                    % QPSK symbols
   sym_map=[1;j;-1;-j];
   Es = pA;
case 8                    % 8PSK symbols
   sym_map=[1;(1+j)/sqrt(2);j;(-1+j)/sqrt(2);-1;(-1-j)/sqrt(2);-j;(1-j)/sqrt(2)];
   Es = pA;
end
 
Eb = Es/K;
na = length(pE);
Ns = 1;              % Number of symbols                   % noise unit variance
 
Es_No = Es/No;
Eb_No = Eb/No;
 
Mt=4;
N=6;

BER = zeros(1,na);
SER = zeros(1,na);
Pseint = zeros(1,na);


for m=1:na
   B=0;
   for loop=1:100000000
       
   bits = round(rand(K,Ns));                        % KxNs matrix of random 0,1 bits
   symbols = graymapPSK(bits);                      % gray code map to symbols
   
%   [h] = manual_channel (6,3,1,PGdb);
   
%Correlation matrix
   R_t = CorrMatrix_interclass(Mt,rhot);
   R_L = CorrMatrix_interclass(L,rhol);

%i.i.d channel
   hw = manual_channel(L,Mt,1,zeros(1,L));
   h = R_t^(1/2)*hw*(R_L.*P)^(1/2);
   
   noise = (randn(1,1) + 1j*randn(1,1))*sqrt(No/2);
   noise1 = (randn(1,1) + 1j*randn(1,1))*sqrt(pE/2);
   
   bit_sequence = sqrt(Es/Ns)*symbols;         
   PSK_seq = bit_sequence(1,:)+1j*bit_sequence(2,:);
   
   %Received signal
     for t=1:Mt
      h_tr(t,:) = conv(h(t,:),conj(fliplr(h(t,:)))/norm(h,'fro'));
      h_tr(t,L) = 0;
     end

   cr = norm(h,'fro')*PSK_seq + noise + noise1;
   
   sd = zeros(2,Ns);
   for n=1:Ns
      [ee ind]=min(abs(sym_map-cr(n)));                    % map to closest symbol
      sd(:,n)= [real(sym_map(ind));imag(sym_map(ind))];    % symbol decision
   end
   %bd=symbols;
   bd = grayunmapPSK(sd,M);
   errors = abs(bd-bits);                % KxNs matrix of bit errors
   symb_err = sign(sum(errors,1));       % one symbol error per column
   SER(m) = sum(symb_err);
   BER(m) = sum(sum(errors))/(K*Ns);
   B = B + SER(m);
   end
   SER(m) = B/(loop*Ns);
end

end






