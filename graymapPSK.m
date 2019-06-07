function symbols = graymapPSK(bits)

K = size(bits,1);
N = size(bits,2);

switch K
case 1             % BPSK
   % maps 0=s1, 1=s0
   % s0 = 1, s1 = -1
   symbols = [bits*2-1;zeros(1,N)];
case 2
   % maps 00=s0, 01=s1, 11=s2, 10=s3
   % s0 = [1;0], s1 = [0;1], s2 = [-1;0], s3 = [0;-1]
case 3
   % maps 000=s0, 001=s1, 011=s2, 010=s3, 110=s4, 111=s5, 101=s6, 100=s7
   % s0 = [1;0],  s1 = 1/sqrt(2)*[1;1],  s2 = [0;1],  s3 = 1/sqrt(2)*[-1;1],
   % s4 = [-1;0], s5 = 1/sqrt(2)*[-1;-1],s6 = [0;-1], s7 = 1/sqrt(2)*[1;-1]
   s = sum(bits,1);
   s_even = [1-bits(1,:)-bits(2,:);bits(2,:)-bits(1,:)];
   s_odd = (1/sqrt(2))*[-1+2*abs(bits(3,:)-bits(1,:));-1+2*abs(bits(3,:)-bits(2,:))];
   symbols = s_even.*([1;1]*(s==0|s==2))+s_odd.*([1;1]*(s==1|s==3));
end