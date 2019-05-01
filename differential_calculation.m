global N wl gama0 B0
gama0 = 10.7083e-4;%C13 parameter
N = 56;
B0 = 401;%magnetic field
wl = 2*pi*gama0*B0;
wh = 1e-3*[30,24,21];
th = pi/180*[133,123,45];
t = (0.05:0.05:10)';
Px = Get_Px(t,wh,th,wl,N);%%这个函数需要传入wl,N才行
M = 2*Px-1;
Px = 
Px_diff = Cal_diff(Px);
plot(t,Px_diff,'color',rand(1,3));
function Px_diff = Cal_diff(Px,t)
%%%算出来的Px_diff为列向量
    switch nargin
        case 1
        Px_diff = zeros(length(Px),1);
        Px_diff(1) = 0;
        for i = 2:length(Px)
            Px_diff(i) = Px(i)-Px(i-1);
        end
        case 2
            Px_diff = zeros(length(Px),1);
            Px_diff(1) = 0;
            for i = 2:length(Px)
                Px_diff(i) = (Px(i)-Px(i-1))/(t(i)-t(i-1));
            end
    end
end
