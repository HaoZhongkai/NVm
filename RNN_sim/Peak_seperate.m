%��peak_seperation3д��һ������,���еĲ����б�����(�󲿷���Ҫ�ں����ڲ���)
%1.����������������:T0,T1
%2.�������ڵĲ���
%3.Ѱ�庯�����ϸ���Ż�
%4.comb���������༰����
%5.Ѱ�����ں��ڸ��������ĵ���
%6.����ٶ���ʱ���ǵȲ�����
%numΪѡ����ĺ˵ĸ���
function wb = Peak_seperate(t,Px,wl,num)
%---------------
%ʹ��findpeaksѰ��
Px1 = Px;
[height,loc,~,~] = findpeaks(Px1,t,'MinPeakProminence',0.1,...
   'Annotate','extents');
%MinPeakDistance',0.005


%---------
%ʹ�ñ任����
%��wl��Χ�ٷ�֮30����
step_size = 0.0001;
T0 = pi/(2*wl)*0.7:step_size:pi/(2*wl)*1.3;
delta = 0.01;
Px2 = zeros(length(loc),1);
for i = 1:length(T0)
    Px2(i) = height'*Discrete_comb(loc,T0(i),delta);
end
%----------------
%���ڴ��ź���ѡȡ����ǿ�ļ���,ɸѡ����Ӧ��������wb
%peak_T�Ƕ�Ӧ������,T_index���ڷ�peak_height�����е��±�
[peak_height,peak_T] = findpeaks(Px2,T0,'MinPeakProminence',1);
[~,T_index] = maxk(peak_height,min(num,length(peak_height)));
T0_signal = peak_T(T_index);
wb = (pi./T0_signal-wl)';


