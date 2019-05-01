%% �ڶ��汾��RNN�Ż�
%��������peak_separate����׼��һ��wb,Ȼ��A,B������һ����Χ�ڽ��е���

%%����׼��
%-------------��������
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.01;
t = (tstep:tstep:tmax)';
%%˥������,���,������,�����źŲ���,���ﲻ�ÿ���˥��
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*20;
%�˵Ĳ���
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
S0 = Kernal(wh,th,wl,N,t);
S0.get_Px();
S0.AddCentralSignal(N_center,wh_center);
P0 = 1-S0.Px;

%% A,B��ֵ�ĳ�ʼ��
param_num = 6;
param_max = [80e-3,60e-3];      %�ֱ��������������С��
num0 = 6;   %��ʼ�����ĺ˵ĸ���
[wb,B_init] = Peak_seperate2(t,wl,N,P0,num0,param_max(1));
%������
% B_init = 0.03*ones(1,num0);
%
init0 = Init_param(wl,wb,B_init,param_num,param_max);
A_init = init0(1:num0);
bound_limit = 0.002;     %Լ��������


%% RNN�ĳ�ʼ���볬��������
RNN_tstep = 0.00001;
% activate_param = [-0.1,0.1];
activate_param = [num0,param_max];
activate_fun = 'SR';
at_fun = Activation_fun2(activate_fun,activate_param);
RNN = RNN2(2*param_num,at_fun);
iterate_num1 = 600;
iterate_num2 = 30;


%% ����RNN�����Ż�
input = init0; 
h = animatedline('color',rand(1,3));
param_exam = zeros(num0,iterate_num1);
output = input;
for i = 1:iterate_num1
%     clearpoints(h);
    for j = 1:iterate_num2
        output  = RNN.Activate(output);
        biase = Constraint_fun(wl,wb,output,num0);
        RNN.Renew(biase);
        output = RNN.feedforward(output,1);
%         addpoints(h,j,output(1));
%         pause(0.001);
    end
    param_exam(:,i) = sqrt((wl+2*pi*output(1:num0)).^2+(2*pi*output(...
        param_num+1:param_num+num0)).^2);
    disp(param_exam(1,i));
    grad = Gradient_loss(P0,output,t,wl,N);
    RNN.Renew(-grad);
%     RNN.Renew(zeros(2*param_num+num0),-grad);
    output = RNN.feedforward(output,RNN_tstep);
    output = RNN.Activate(output);
    loss = Get_loss(P0,output(1:2*param_num),t,wl,N);
    addpoints(h,i,loss);
%     addpoints(g,i,output(3));
    pause(0.0001);
end

%% ������������
A1 = output(1:param_num);
B1 = output(param_num+1:param_num*2);
P1 = (1-Get_M2(t,A1,B1,wl,N))/2;
figure
hold on;
scatter(init0(1:param_num),init0(param_num+1:2*param_num));
scatter(A0,B0,'MarkerEdgeColor','r','SizeData',60);
title('��ʼ��A,B�ķֲ�')
figure
hold on;
grid on;
scatter(output(1:param_num),output(param_num+1:2*param_num),'SizeData',30);
scatter(A0,B0,'MarkerEdgeColor','r','SizeData',60);
scatter(init0(1:param_num),init0(param_num+1:2*param_num),'MarkerEdgeColor','g');
title('ĩ̬��A,B�ķֲ�')
figure
hold on;
plot(t,P0,t,P1);


%     T_conf = T_confMatrix(wl,input,num0);
%     biase = zeros(2*param_num+num0,1);
%     RNN.Renew(T_conf,biase);
%         E_conf = -output'*T_confMatrix(wl,output,num0)*output;





