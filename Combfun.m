%�������Դ��м��ĺ����ų�һ�е����ź�
%T0Ϊ��ʼ����,delta����ΪGauss��Ŀ��
%����yΪ������
function y = Combfun(t,T0,delta)
N = floor((t(end)/T0+1)/2);%�������
y = zeros(length(t),1);
for i = 1:N
    y = y+exp(-((t-(2*i-1)*T0)/delta).^2);%��˹����
%     y = y+1./(((t-(2*i-1)*T0)/delta).^2+1);
end
