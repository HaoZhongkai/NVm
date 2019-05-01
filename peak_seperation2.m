%-------------
%另外一种基于原数据对峰进行尝试拟合的方法
%使用一个变换将原函数有峰的地方挑选出来
%输入数据为已经去掉基线的数据
%核的基本参数
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
%衰减周期,误差,脉冲数,中心信号参数,这里不用考虑衰减
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*15;
%核的参数
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
% wh = 1e-3*[33];
% th = pi/180*[132];
%计算核的参数
S = Kernal(wh,th,wl,N,t);
S.get_Px();
S.AddCentralSignal(N_center,wh_center);
% S.Modulate(Ta,Tb);
% S.Addnoise(e);
Px = S.Px;
figure 
axis([0 tmax 0 1])
plot(t,Px);


%---------------------
%考虑用一个函数来对f进行变换
%尝试使用几个Gauss函数的核
T0 = 0.4:0.001:0.7;
delta = 0.0001;
Px1 = zeros(length(T0),1);
for i = 1:length(T0)
    Peakfun = Combfun(t,T0(i),delta);
    Px1(i) = (1-Px)'*Peakfun;
end
figure
plot(T0,Px1);

%---------
%下一步，将所有local_peak的地方



