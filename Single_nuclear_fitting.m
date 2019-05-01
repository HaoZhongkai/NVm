%%%----------------------------
%���һ��:�������źŵ����
%������ʹ�ø��ַ����Ը��ֲ����ĺ˽������
%���������ѵĵط�������ο˷����ܴ��ںܴ����ĺ˵����
%�˵Ļ�������
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
%˥������,���,������,�����źŲ���,���ﲻ�ÿ���˥��
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*20;
%�˵Ĳ���
% wh = 1e-3*[83.8,47,55,19,33,25.1];
% th = pi/180*[21,30,54,133,132,51];
wh = 1e-3*[83];
th = pi/180*[21];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
%����˵Ĳ���
S = Kernal(wh,th,wl,N,t);
S.get_Px();
% S.AddCentralSignal(N_center,wh_center);
% S.Addnoise(e);
Px = S.Px;
[T_signal,t_signal,Px_signal] = Peak_seperate(t,Px);

%--------------
%�����ź�
figure 
hold on;
plot(t,Px);
axis([0 tmax 0 1]);

%�Ժ��ٴ������ź�ʱ����signal��Ϊ1-signal
Px_signal = 1-Px_signal;



%-------------------
%���濼��������ĳ������Ϊĳ��ԭ��ȱʧ�����Ǳ�Ϊ�����ֵ��ʱ��
%�ȿ���һ�������߶�ʧ�����
index_lost = randi(length(t_signal));
% Px_signal(index_lost) = 0;
%��������õ��ķ��ͼ��
scatter(t_signal,1-Px_signal,'MarkerEdgeColor','r');



%---����A�Ľ���ֵ��beta=wb/wl��׼ȷֵ,������ԭ���ĺ˵�beta׼ȷֵ
Ac = 1/(2*sum(t_signal)/length(t_signal)^2)-wl/pi;
beta0 = S.get_wb./wl;
beta = pi/(sum(t_signal)*wl/length(t_signal)^2)-1;
A_t = (2*(1:length(t_signal))'-1)*pi./t_signal-2*wl;

%--------------����1:���ȡ��
%����B1,�����ź�,��Bɨ��һ��,������ʧ����
%ÿ�����ź���ȥ��������������

% len = length(t_signal);
% B_best = zeros(len*(len-1)/2,1);
% iterate = 1;
% 
% for loop1 = 1:length(t_signal)
%     for loop2 = loop1+1:length(t_signal)
% %         Px_signal_copy = Px_signal;%��ʱ����������copy
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

%----------------------����2:������Ϸ�
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


