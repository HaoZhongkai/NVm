B0 = 400;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
N = 100;
wh = 2*pi*83.8*1e-4;
th = 21*pi/180;
A = wh*cos(th);
B = wh*sin(th);
wl = gama0*B0;
wb = sqrt((A+wl)^2+B^2);
mz = (A+wl)/wb;
mx = B/wb;
tmax = 20;
t = 0:0.1:tmax;
M = ones(1,length(t));
kmax = floor((tmax*(2*wl+A)/pi-1)/2);
k = 1:kmax;
tao = (2*k+1)*pi/(2*wl+A);
dealt = -0.2:0.005:0.2;
M0 = zeros(1,length(dealt)*kmax);
tao0 = zeros(1,length(dealt)*kmax);
for i = k
    dk = (2*i+1)*pi*dealt;
    tao0(1,(i-1)*length(dealt)+1:i*length(dealt)) = dk+tao(i);
    M0(1,(i-1)*length(dealt)+1:i*length(dealt)) = ...
        1-2./(1+dk.^2./mx.^2).*sin(N/2*sqrt(mx^2+dk.^2)).^2;
end
tao0 = reshape(tao0,[1,kmax*length(dealt)]);
M0 = reshape(M0,[1,kmax*length(dealt)]);
Px0 = (1+M0)/2;
Px = (1+M)/2;
figure
hold on;
plot(t,Px,'color',rand(1,3));
plot(tao0,Px0,'color',rand(1,3));

