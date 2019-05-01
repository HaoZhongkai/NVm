%在一有界区域内进行最小搜索
%尝试使用随机方法进行优化
%因为该函数并非最小二乘函数，应使用凸优化方法来优化
function [Ta,Tb] = LMfun(Ta_max,Tb_max,t,BaseIndex,Px_data,N)
%将网格分为5等分随机初始化
N_grid = 5;
Ta0 = zeros(1,N_grid);
Tb0 = zeros(1,N_grid);
for i = 1:N_grid
    Ta0(i) = Ta_max/N_grid*rand+(i-1)*Ta_max/N_grid;
    Tb0(i) = Tb_max/N_grid*rand+(i-1)*Tb_max/N_grid;
end

%--------
%初始化误差，使用LM方法最小化
%LM方法跳出的原则如下:
%   1.误差迭代几乎不再减小
%   2.超过指定循环次数
S = 0;
S_best = 1e10;      %S初值设为很大
gama = 1;
for i = 1:N
    for j = 1:N
        J = Jacobifun(Ta0(i),Tb(j),Px_data,t,BaseIndex,N);
       
        
        
        
