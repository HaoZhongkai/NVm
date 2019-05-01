B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
tmax = 10;
tstep = 0.001;
t = (0.01:tstep:tmax)';
Ta = 5e3;
Tb = 5e3;
%set parameter...
wh = 1e-3*[83.8,47];
th = pi/180*[21,30];
% wh = 1e-3*[30];
% th = pi/180*[133];
% wh_c = 1e-3*20*rand(1,60);
% th_c = pi*rand(1,60);
% wh = [wh0,wh_c];
% th = [th0,th_c];
N = 56; 
wl = 2*pi*gama0*B0;
%set parameter...
A = wh.*cos(th);
B = wh.*sin(th);
Px = Get_Px(t,wh,th,wl,N);
%任意设置一个数据点，判断它是不是信号
wb = sqrt((2*pi*A+wl).^2+(2*pi*B).^2);
T0 = pi/(wl+wb(2));
%在周N00个点寻找局部峰
N00 = 10;
% i0 = randi([40, 70],1);
i1 = 0;
while(floor((2*i1+1)*T0*1e3)<length(t))
    i1 = i1+1;
end
i1 = i1-1;
signal_index = floor(T0*1e3*(2*(0:i1)+1));
[signal,signal_index] = local_peak(Px,signal_index,N00);
figure
hold on;
plot(1:length(t),Px,'color','b');
scatter(signal_index,signal,30,'filled');
% figure
% hold on;
% scatter(1:length(signal_index),signal_index);




