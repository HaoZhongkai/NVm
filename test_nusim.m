%--------------
%调试的时候看Px的值
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
Ta = 5e3;
Tb = 5e3;
%set parameter...
% wh0 = 1e-3*[83.8,47,55,19,33,25.1];
% th0 = pi/180*[21,30,54,133,132,51];
wh0 = 1e-3*[50,40,32];
th0 = pi/180*[32,46,19];
% wh_c = 1e-3*20*rand(1,60);
% th_c = pi*rand(1,60);
wh_c = [];
th_c = [];
wh = [wh0,wh_c];
th = [th0,th_c];
A = wh.*cos(th);
B = wh.*sin(th);
N = 56; 
wl = 2*pi*gama0*B0;
beta = t*wl;
M = ones(length(t),1);
ep = 0.01*randn(length(t),1);
for i = 1:length(A)
    wb = sqrt((2*pi*A(i)+wl)^2+(2*pi*B(i))^2);
    mz = (2*pi*A(i)+wl)/wb;
    mx = 2*pi*B(i)/wb;
    alpha = t*wb;
    Cphi = cos(alpha).*cos(beta)-mz*sin(alpha).*sin(beta);
    phi = acos(Cphi);
    n01 = mx^2*(1-cos(alpha)).*(1-cos(beta))./(1+Cphi);
    M = M.*(1-n01.*sin(N*phi/2).^2);
end
Px = (1+M)/2;
% Px = 1/2*M.*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb)+ep;
for i = 1:length(t)
    if Px(i)>1
         Px(i) = Px(i)-2*ep(i);
    end
end
figure
hold on;
axis([0 tmax 0 1]);
plot(t,Px,'color',rand(1,3));
