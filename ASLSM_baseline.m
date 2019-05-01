%% 使用非对称最小二乘法来求基线，目标是使用非对称权重将基线附近的点数最大化
% t为输入的时间,base_index为属于基线的点的下标,ep为噪声的振幅
% 求解的方法为梯度下降(上升)
%weight为带
function coef = ASLSM_baseline(t,Px,fun,Base_index,ep,coef0,weight,rate)
    t_base = t(Base_index);
    Px_base = Px(Base_index);
    coef = coef0;
%     h = animatedline('color',rand(1,3));
    i = 0;
    pointnum = 1;
    pointnum0 = 0;
    while(abs(pointnum-pointnum0)/pointnum>0.000&&i<4000)
        Px_fit = fun(coef,t_base);
        pointnum0 = pointnum;   %上次的值
        pointnum = Lossfun(Px_base,Px_fit,weight,ep);
        grad = Grad_num(fun,Px_base,t_base,weight,coef,ep);
        coef = coef+rate*grad;
%         addpoints(h,i,pointnum);
%         pause(0.001);
        i = i+1;
    end
    %调试部分
end
%连续化的点的个数函数
function pointnum = Lossfun(Px_base,Px_fit,weight,ep)
    pointnum = sum((Px_base>=Px_fit).*exp((Px_fit-Px_base)/ep)+weight*(Px_base...
        <Px_fit).*exp((Px_base-Px_fit)/ep));
end
%计算梯度grad为行向量(为了匹配coef少有的情况)
function grad = Grad_num(fun,Px_base,t_base,weight,coef,ep)
    step = 1e-6;%固定参数:求导的步长
    Px_fit0 = fun(coef,t_base);
    pointnum0 = Lossfun(Px_base,Px_fit0,weight,ep);
    coef1 = coef+[step,0];
    coef2 = coef+[0,step];
    Px_fit1 = fun(coef1,t_base);
    Px_fit2 = fun(coef2,t_base);
    pointnum1 = Lossfun(Px_base,Px_fit1,weight,ep);
    pointnum2 = Lossfun(Px_base,Px_fit2,weight,ep);
    grad = [pointnum1-pointnum0,pointnum2-pointnum0]/ep;
end



