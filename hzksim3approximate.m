mx = 0.04;
N = 90;
dx = -0.5:0.001:0.5;
M = 1-2./(1+dx.^2/mx^2).*sin(N/2*sqrt(mx^2+dx.^2)).^2;
P = (1+M)/2;
plot(dx,P);