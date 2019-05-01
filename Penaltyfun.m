%����С���˷�ȷ����������Сֵ
%���Ϊ����,���ڸ���������,����Ҫƽ����ÿ����
%pΪȨ��,������������·�ʱȨ�ؽ�С,����������Ϸ�ʱ���ϴ�
function S = Penaltyfun(Ta,Tb,Px_data,t,BaseIndex,N)
p = 0.05;
Ta = Ta*1e3;
Tb = Tb*1e3;
Px0 = 1/2*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb);
S = 0;
for i = BaseIndex
    S = S+p*max(Px0(i)-Px_data(i),0)^2+(1-p)*min(Px0(i)-Px_data(i),0)^2;
end