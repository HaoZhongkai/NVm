%-------------产生数据
%--------------------------------------------------------------------------
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
%%衰减周期,误差,脉冲数,中心信号参数,这里不用考虑衰减
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*20;
%核的参数
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
% wh = 1e-3*[33];
% th = pi/180*[132];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
%计算核数据
S = Kernal(wh(1),th(1),wl,N,t);
S.get_Px();
Px0 = 1-S.Px;
S.AddCentralSignal(N_center,wh_center);
% S.Addnoise(e);
Px = 1-S.Px;
S1 = Kernal(wh(1),th(2),wl,N,t);
S1.get_Px();
Px1 = 1-S1.Px;
r0 = Px0'*Px
r1 = Px1'*Px

