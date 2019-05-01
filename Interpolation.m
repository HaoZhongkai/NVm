function y0 =Interpolation(y, t, t0)
n = length(t);
for i = 1:n-1
    if t0>=t(i)&&t0<t(i+1)
        y0 = y(i)+(y(i+1)-y(i))*(t0-t(i))/(t(i+1)-t(i));
    end
end
        