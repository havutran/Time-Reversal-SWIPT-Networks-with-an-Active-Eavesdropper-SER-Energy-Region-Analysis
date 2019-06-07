function bits = grayunmapPSK(symbols,M)

K = log2(M);
N = size(symbols,2);

switch K
case 1             % BPSK
   % maps 0=s1, 1=s0
   % s0 = 1, s1 = -1
   bits = 0.5*(symbols(1,:)+1);
   case 2
   % maps 00=s0, 01=s1, 11=s2, 10=s3
   % s0 = [1;0], s1 = [0;1], s2 = [-1;0], s3 = [0;-1]
   bits = [0.5*(1-symbols(1,:)-symbols(2,:));0.5*(1+symbols(2,:)-symbols(1,:))];
   case 3
   % maps 000=s0, 001=s1, 011=s2, 010=s3, 110=s4, 111=s5, 101=s6, 100=s7
   % s0 = [1;0],  s1 = 1/sqrt(2)*[1;1],  s2 = [0;1],  s3 = 1/sqrt(2)*[-1;1],
   % s4 = [-1;0], s5 = 1/sqrt(2)*[-1;-1],s6 = [0;-1], s7 = 1/sqrt(2)*[1;-1]
   bits_even = [0.5*(1-symbols(1,:)-symbols(2,:));0.5*(1+symbols(2,:)-symbols(1,:));...
         abs(symbols(2,:))];
   bits_odd = 0.5*[1-symbols(2,:)*sqrt(2);1-symbols(1,:)*sqrt(2);sqrt(2)*abs(symbols(1,:)+symbols(2,:))];
   s = abs(symbols(1,:))>0.7 & abs(symbols(1,:))<0.8;
   bits = bits_even.*([1;1;1]*(~s))+bits_odd.*([1;1;1]*(s));
end