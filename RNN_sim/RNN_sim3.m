%% 第三版本的RNN优化
%首先利用peak_separate函数准备一组wb,然后将A,B限制在一个范围内进行迭代

%%数据准备
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

%% A,B的值的初始化
param_num = 12;
param_max = [90e-3,50e-3];      %分别用来整定主峰和小峰
num0 = 6;   %初始整定的核的个数
wb = Peak_seperate(t,P0,wl,6);
init0 = Init_param(wl,wb,param_num,param_max);
bound_limit = 0.002;     %约束条件的


%% RNN的初始化与超参数整定
RNN_tstep = 0.0001;
activate_param = [-0.1,0.1];
activate_fun = 'SR';
at_fun = Activation_fun(activate_fun,activate_param);
RNN = RNN2(2*param_num+2*num0,at_fun);
iterate_num1 = 600;
iterate_num2 = 400;


%% 利用RNN进行优化
input = [init0;zeros(2*num0,1)]; %补上的为松弛变量
h = animatedline('color',rand(1,3));
g = animatedline('color',rand(1,3));
output = input;
for i = 1:iterate_num1
    for j = 1:iterate_num2
%     T_conf = T_confMatrix(wl,input,num0);
%     biase = zeros(2*param_num+num0,1);
        biase = Constraint_fun(wl,wb,output,num0,bound_limit);
        RNN.Renew(biase);
%     RNN.Renew(T_conf,biase);
        output = RNN.feedforward(output,1);
        output  = RNN.Activate(output);
%         E_conf = -output'*T_confMatrix(wl,output,num0)*output;
%         addpoints(h,j,E_conf)
    end
    grad = [Gradient_loss2(P0,output,t,wl,N,param_num);zeros(num0,1)];
    RNN.Renew(-grad);
%     RNN.Renew(zeros(2*param_num+num0),-grad);
    output = RNN.feedforward(output,RNN_tstep);
    output = RNN.Activate(output);
    loss = Get_loss(P0,output(1:2*param_num),t,wl,N);
    addpoints(h,i,loss);
%     addpoints(g,i,output(3));
    pause(0.001);
end

%% 输出结果并分析
A1 = output(1:param_num);
B1 = output(param_num+1:2*param_num);
P1 = (1-Get_M2(t,A1,B1,wl,N))/2;
figure
hold on;
grid on;
scatter(output(1:param_num),output(param_num+1:2*param_num),'SizeData',30);
scatter(A0,B0,'MarkerEdgeColor','r','SizeData',40);
figure
hold on;
plot(t,P0,t,P1);




