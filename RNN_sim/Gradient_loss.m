%计算给定的S和目标P0在L2范数下对参数A,B的梯度函数
%输入是一个标准信号和一个拟合用的kernal
%输出output为一个列向量
function grad = Gradient_loss(P0,param,t,wl,N)
    ep = 1e-7;
    tstep = t(2)-t(1);
    M0 = 1-2*P0;
    len = length(param)/2;
    A = param(1:len);
    B = param(len+1:end);
    M1 = Get_M2(t,A,B,wl,N);
    grad = zeros(2*len,1);
    biase = (M1-M0)'*(M1-M0);
    for i = 1:len
        A2 = A;
        B2 = B;
        A2(i) = A(i)+ep;
        B2(i) = B(i)+ep;
        M21 = Get_M2(t,A2,B,wl,N);
        M31 = Get_M2(t,A,B2,wl,N);
        grad(i) = tstep*((M21-M0)'*(M21-M0)-biase)/(4*ep);
        grad(i+len) = tstep*((M31-M0)'*(M31-M0)-biase)/(4*ep);
        A2(i) = A(i)-ep;
        B2(i) = B(i)-ep;
        M22 = Get_M2(t,A2,B,wl,N);
        M32 = Get_M2(t,A,B2,wl,N);
        grad(i) = tstep*((M21-M0)'*(M21-M0)-(M22-M0)'*(M22-M0))/(8*ep);
        grad(i+len) = tstep*((M31-M0)'*(M31-M0)-(M32-M0)'*(M32-M0))/(8*ep);
    end
end