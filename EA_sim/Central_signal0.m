%本程序用于定性观察

B0 = 675;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
wh = [];
th = [];
N = 24;
% t = (0.01:0.01:10)';
num0 = 100;
wh_c = 1e-3*14*rand(1,num0);
th_c = pi*rand(1,num0);
wh_in = [wh,wh_c];
th_in = [th,th_c];
Px0 = Get_Px(t,wh_in,th_in,wl,N);
plot(t,Px,t,Px0);
axis([0 t(end) 0 1]);