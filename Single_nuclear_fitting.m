%%%----------------------------
%最后一步:单个核信号的拟合
%将尝试使用各种方法对各种参数的核进行拟合
%这里最困难的地方在于如何克服可能存在很大误差的核的拟合
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
wh_center = 1e-3*20;
%核的参数
% wh = 1e-3*[83.8,47,55,19,33,25.1];
% th = pi/180*[21,30,54,133,132,51];
wh = 1e-3*[83];
th = pi/180*[21];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
%计算核的参数
S = Kernal(wh,th,wl,N,t);
S.get_Px();
% S.AddCentralSignal(N_center,wh_center);
% S.Addnoise(e);
Px = S.Px;
[T_signal,t_signal,Px_signal] = Peak_seperate(t,Px);

%--------------
%单核信号
figure 
hold on;
plot(t,Px);
axis([0 tmax 0 1]);

%以后再处理单核信号时都将signal变为1-signal
Px_signal = 1-Px_signal;



%-------------------
%下面考虑谱线中某个峰因为某种原因缺失，或是变为错误的值的时候
%先考虑一条峰谱线丢失的情况
index_lost = randi(length(t_signal));
% Px_signal(index_lost) = 0;
%补画分离得到的峰的图像
scatter(t_signal,1-Px_signal,'MarkerEdgeColor','r');



%---计算A的近似值即beta=wb/wl的准确值,并计算原来的核的beta准确值
Ac = 1/(2*sum(t_signal)/length(t_signal)^2)-wl/pi;
beta0 = S.get_wb./wl;
beta = pi/(sum(t_signal)*wl/length(t_signal)^2)-1;
A_t = (2*(1:length(t_signal))'-1)*pi./t_signal-2*wl;

%--------------方法1:随机取点
%给定B1,计算信号,将B扫描一遍,计算损失函数
%每次在信号中去除两个点进行拟合

% len = length(t_signal);
% B_best = zeros(len*(len-1)/2,1);
% iterate = 1;
% 
% for loop1 = 1:length(t_signal)
%     for loop2 = loop1+1:length(t_signal)
% %         Px_signal_copy = Px_signal;%临时用来操作的copy
% %         t_signal_copy = t_signal;
%         Px_signal_copy = Px_signal([loop1,loop2]);        
%         t_signal_copy = t_signal([loop1,loop2]);
%         B_arr = linspace(0,0.05,100);
%         Biase = zeros(length(B_arr),1);
%         for i = 1:length(B_arr)
%             signal1 = Single_peakdepth(B_arr(i),wl,beta0,N,t_signal_copy);
%             Biase(i) = (Px_signal_copy-signal1)'*(Px_signal_copy-signal1);
%         end
%         [Biase_min,index_mloss] = min(Biase);
%         B_best(iterate) = B_arr(index_mloss);
%         iterate = iterate+1;
%     end
% end
% figure
% hold on;
% scatter(1:length(B_best),B_best);
%--------------------------------

%----------------------方法2:迭代拟合法
B_arr = linspace(0,0.05,100);
Biase = zeros(length(B_arr),1);
for i = 1:length(B_arr)
    signal1 = Single_peakdepth(B_arr(i),wl,beta0,N,t_signal);
    Biase(i) = (Px_signal-signal1)'*(Px_signal-signal1);
end
[Biase_min,index_mloss] = min(Biase);
B_best = B_arr(index_mloss)
signal1 = Single_peakdepth(B_best,wl,beta,N,t_signal);
scatter(t_signal,1-signal1);

% figure
% hold on;
% scatter(t_signal,Px_signal,'MarkerFacecolor','r');
% scatter(t_signal,signal1);
% axis([0 tmax 0 1]);


