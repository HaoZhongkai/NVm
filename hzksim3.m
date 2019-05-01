%根据A,B观察理论图像
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
tmax = 100;
tstep = 0.001;
t = (0.01:tstep:tmax)';
Ta = 5e3;
Tb = 5e3;
N = 56; 
A = 1e-3*[20,10];
B = 1e-3*[10,10];
%set parameter...
% wh0 = 1e-3*[83.8,47,55,19,33,25.1];
% th0 = pi/180*[21,30,54,133,132,51];
% wh_c = 1e-3*20*rand(1,60);
% th_c = pi*rand(1,60);
% wh = [wh0,wh_c];
% th = [th0,th_c];
%set parameter...
% A = wh.*cos(th);
% B = wh.*sin(th);
wl = 2*pi*gama0*B0;
Px1 = Get_Px2(t,A(1),B(1),wl,N);
Px2 = Get_Px2(t,A(2),B(2),wl,N);
%plot Px
% Px = 1/2*M+1/2;
% Px = 1/2*M.*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb)+ep;
figure
hold on;
axis([0 tmax 0 1]);
plot(t,Px1,'color','r');
plot(t,Px2,'color','b');