%用最小二乘法确定罚函数最小值
%输出为误差函数,对于给定的数据,不需要平均到每个点
%p为权重,点在拟合曲线下方时权重较小,在拟合曲线上方时误差较大
function S = Penaltyfun(Ta,Tb,Px_data,t,BaseIndex,N)
p = 0.05;
Ta = Ta*1e3;
Tb = Tb*1e3;
Px0 = 1/2*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb);
S = 0;
for i = BaseIndex
    S = S+p*max(Px0(i)-Px_data(i),0)^2+(1-p)*min(Px0(i)-Px_data(i),0)^2;
end