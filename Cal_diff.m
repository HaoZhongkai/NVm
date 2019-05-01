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