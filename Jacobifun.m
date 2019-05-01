%计算LM方法中的Jacobi矩阵
function J = Jacobifun(Ta,Tb,Px_data,t,BaseIndex,N)
step = 1e-3;
J = zeros(length(BaseIndex),2);
J(:,1) = (Penaltyfun(Ta+step,Tb,Px_data,t,BaseIndex,N)-Penaltyfun...
    (Ta,Tb,Px_data,t,BaseIndex,N))/step;
J(:,2) = (Penaltyfun(Ta,Tb+step,Px_data,t,BaseIndex,N)-Penalyfun...
    (Ta,Tb,Px_data,t,BaseIndex,N))/step;
end
