%��һ�н������ڽ�����С����
%����ʹ��������������Ż�
%��Ϊ�ú���������С���˺�����Ӧʹ��͹�Ż��������Ż�
function [Ta,Tb] = LMfun(Ta_max,Tb_max,t,BaseIndex,Px_data,N)
%�������Ϊ5�ȷ������ʼ��
N_grid = 5;
Ta0 = zeros(1,N_grid);
Tb0 = zeros(1,N_grid);
for i = 1:N_grid
    Ta0(i) = Ta_max/N_grid*rand+(i-1)*Ta_max/N_grid;
    Tb0(i) = Tb_max/N_grid*rand+(i-1)*Tb_max/N_grid;
end

%--------
%��ʼ����ʹ��LM������С��
%LM����������ԭ������:
%   1.�������������ټ�С
%   2.����ָ��ѭ������
S = 0;
S_best = 1e10;      %S��ֵ��Ϊ�ܴ�
gama = 1;
for i = 1:N
    for j = 1:N
        J = Jacobifun(Ta0(i),Tb(j),Px_data,t,BaseIndex,N);
       
        
        
        
