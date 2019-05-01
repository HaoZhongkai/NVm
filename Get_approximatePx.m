function [t_appro,appro_Px] = Get_approximatePx(wl,A,n,mx,N,delta,N0)
%输入n为需要画的n个周期的近似峰的个数,N0为每个峰附近点数,delta为画的峰的范围
t_appro = zeros(N0*n,1);
delta = linspace(-delta,delta,N0);
appro_Px = zeros(N0*n,1);
for i = 1:n
    tao_k = (2*i-1)*pi/(2*wl+A);
    delta_k = (2*i-1)*pi*delta;
    t_appro((i-1)*N0+1:i*N0) = tao_k*(1+delta);
    appro_Px((i-1)*N0+1:i*N0) = 1-sin(N*sqrt(mx^2+delta_k.^2)/2).^2./(...
        1+delta_k.^2/(mx.^2));
    appro_Px((i-1)*N0+1) = 1;
    appro_Px(i*N0) = 1;
end

