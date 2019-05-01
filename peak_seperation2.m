%-------------
%����һ�ֻ���ԭ���ݶԷ���г�����ϵķ���
%ʹ��һ���任��ԭ�����з�ĵط���ѡ����
%��������Ϊ�Ѿ�ȥ�����ߵ�����
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
wh_center = 1e-3*15;
%�˵Ĳ���
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
% wh = 1e-3*[33];
% th = pi/180*[132];
%����˵Ĳ���
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
%������һ����������f���б任
%����ʹ�ü���Gauss�����ĺ�
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
%��һ����������local_peak�ĵط�



