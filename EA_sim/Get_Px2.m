function Px = Get_Px2(t,A,B,wl,N)
%这个函数要求输入A和B
beta = t*wl;
M = ones(length(t),1);
for i = 1:length(A)
    wb = sqrt((2*pi*A(i)+wl)^2+(2*pi*B(i))^2);
    mz = (2*pi*A(i)+wl)/wb;
    mx = 2*pi*B(i)/wb;
    alpha = t*wb;
    Cphi = cos(alpha).*cos(beta)-mz*sin(alpha).*sin(beta);
    phi = acos(Cphi);
    n01 = mx^2*(1-cos(alpha)).*(1-cos(beta))./(1+Cphi);
    M = M.*(1-n01.*sin(N*phi/2).^2);
end
Px = 1/2*M+1/2;