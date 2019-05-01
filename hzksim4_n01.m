wl = 2.698;
A = 1e-3*20;
B = 1e-3*25;
th = pi/180*75;
gama = 0.029;
step = 0.0005;
alpha_max = 500;
t = 0:step:alpha_max;
alpha = wl*t;
N = 89;
c = sqrt(gama^2+2*gama*cos(th)+1);
beta = c*alpha;
mx = gama*sin(th)/sqrt(gama^2+2*gama*cos(th)+1);
mz = (gama*cos(th)+1)/sqrt(gama^2+2*gama*cos(th)+1);
n01 = mx^2*(1-cos(alpha)).*(1-cos(beta))./(1+cos(alpha).*cos(beta)...
    -mz*sin(alpha).*sin(beta));
num_period = 2000;
alpha1 = (1:2:2*num_period-1)*pi/(1+c);
beta1 = c*alpha1;
Cphi = cos(alpha+beta)+(1-mz)/2*(cos((c-1)*alpha)-cos((c+1)*alpha));
phi = acos(Cphi);
phi1 = acos((1-mz)*cos((c-1)*alpha/2).^2+cos(alpha+beta));
%%%phi2,phi1为近似后的phi
phi2 = pi-abs(gama*sin(th)*cos((c-1)*alpha/2));
y1 = sin(N*(phi)/2).^2;
y2 = sin(N*(phi1)/2).^2;
M1 = n01.*y1;
M2 = n01.*y2;
figure;
hold on;
plot(alpha,M1);
% plot(alpha,M2,'color','r');
% T = 2*pi/(gama*cos(th))
% plot(alpha1,y1*2);
% plot(alpha,n01,'color','b');

