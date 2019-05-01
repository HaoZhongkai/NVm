%最简单的loss fun，方便主函数调试与评估性能用
function loss = Get_loss(P0,param,t,wl,N)
    M1 = Get_M2(t,param(1:end/2),param(end/2+1:end),wl,N);
    P1 = (1-M1)/2;
    loss = (P1-P0)'*(P1-P0)*(t(2)-t(1)); 
end