%% 使用RNN进行参数优化
%% 数据准备
%-------------产生数据
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.01;
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
A0 = wh.*cos(th);
B0 = wh.*sin(th);
S0 = Kernal(wh,th,wl,N,t);
S0.get_Px();
S0.AddCentralSignal(N_center,wh_center);
P0 = 1-S0.Px;

%% RNN参数设计与初始化
RNN_tstep = 0.00002;
b_max = 20e-3;
param_num = 6;
activate_param = [-0.1,0.1];
activate_fun = 'SR';
at_fun = Activation_fun1(activate_fun,activate_param);
RNN = RNN1(param_num,at_fun);
iterate_num = 600;
% Init_param = b_max*[(2*rand(param_num,1)-1);rand(param_num,1)];
Init_param = [A0';b_max*(2*rand(param_num-6,1)-1);B0'+1e-2;...
    b_max*rand(param_num-6,1)];
% Init_param = [A0';B0'];

%% 利用RNN进行参数优化
h = animatedline('color',rand(1,3));
input = Init_param;
for i = 1:iterate_num
    grad = Gradient_loss(P0,input,t,wl,N);
    RNN.Renew(grad);
    output = RNN.feedforward(input,RNN_tstep);
    output = RNN.Activate(output);
    input = output;
    loss = Get_loss(P0,output(1:2*param_num),t,wl,N);
    %绘图
    addpoints(h,i,loss);
    pause(0.001);
end
%% 输出结果并分析
A1 = output(1:param_num);
B1 = output(param_num+1:end);
P1 = (1-Get_M2(t,A1,B1,wl,N))/2;
figure
hold on;
grid on;
scatter(output(1:param_num),output(param_num+1:end),'SizeData',30);
scatter(A0,B0,'MarkerEdgeColor','r','SizeData',40);
figure
hold on;
plot(t,P0,t,P1);

