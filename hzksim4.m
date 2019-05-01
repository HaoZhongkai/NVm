%观察单个核的调制信号sin
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
A = 1e-3*50;
B = 1e-3*50;
tmax = 10;
tstep = 0.001;
t = (0.01:tstep:tmax)';
M = ones(length(t),1);
Ta = 5e3;
Tb = 5e3;
N = 89;
wl = 2*pi*gama0*B0;
beta = t*wl;
wb = sqrt((2*pi*A+wl)^2+(2*pi*B)^2);
alpha = t*wb;
mz = (2*pi*A+wl)/wb;
mx = (2*pi*B)/wb;
Cphi = cos(alpha).*cos(beta)-mz*sin(alpha).*sin(beta);
phi = acos(Cphi);
Q = sin(N*phi/2).^2;
n01 = mx^2*(1-cos(alpha)).*(1-cos(beta))./(1+Cphi);
Px = 1-n01.*Q/2;
figure
plot(t,Px,'color','r');
axis([0 tmax 0 1]);
figure
hold on;
% plot(t,cos(alpha),'g');
% plot(t,cos(beta),'r');
% plot(t,Q,'color','b');
plot(t,n01,'color','b');
% s = scatter(t,Q);
% plot(t,Q);
% s.SizeData = 10;
