%�����Ӧ��B(������),wl,beta(��wb/wl),N,t,�����Ӧ��B�ķ�ֵ��С(1-Px)(������)
function signal = Single_peakdepth(B,wl,beta,N,t_signal)

%���ϼн���ĵ���
alpha = wl*beta.*t_signal;
beta_t = wl*t_signal;%��������������һ��beta
mx = 2*pi*B./(beta*wl);
mz = sqrt(1-mx^2);
n01 = mx^2*(1+cos(alpha)).*(1+cos(beta_t))./(1+cos(alpha).*cos(beta_t)-...
    mz*sin(alpha).*sin(beta_t));
phi = acos(cos(beta*t_signal*wl).*cos(t_signal*wl)-...
    sqrt(1-(2*pi*B/(beta*wl)).^2).*sin(wl*beta*t_signal).*sin(wl*t_signal));
% phi = acos(-1+(1-sqrt(1-(2*pi*B./(beta*wl)).^2))*(cos((1-beta)...
%     *t_signal/2).^2).^2);
signal = n01.*sin(N*phi/2).^2/2;
end