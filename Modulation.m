function y = Modulation(Px,t,N,Ta,Tb)
y = 1/2*(2*Px-1).*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb);

