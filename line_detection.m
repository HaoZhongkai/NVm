B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
Ta = 5e3;
Tb = 5e3;
%set parameter...
wh = 1e-3*[83.8,47];
th = pi/180*[21,30];
N = 56; 
wl = 2*pi*gama0*B0;
wb = sqrt((2*pi*wh).^2+2*2*pi*wl*wh.*cos(th)+wl^2);
A = wh.*cos(th);
B = wh.*sin(th);
Px = Get_Px(t,wh,th,wl,N);
Px2 = 1-Px;
%-----------------------
%使用matlab寻峰算法寻峰
[height,loc,pk_w,pk_p] = findpeaks(Px2,t,'MinPeakProminence',0.03,...
    'MinPeakDistance',0.05,'Annotate','extents');
[k_line,b_line] = line_detect(loc);
figure
hold on;
scatter(k_line,b_line);      
figure
hold on;
scatter(loc,1-height,'MarkerEdgeColor','r');
plot(t,Px);



%-----------------------------
%line detection
%任意取出两个点，记录其直线然后按斜率截距画图
function [k_line,b_line] = line_detect(loc)
    n = length(loc);
    k_line = zeros(1,n*(n-1)/2);
    b_line = zeros(1,n*(n-1)/2);
    loop = 1;
    for i = 1:n
        for j = i+1:n
            k_line(loop) = loc(j)-loc(i);
            b_line(loop) = loc(i);
            loop = loop+1;
        end
    end
end
