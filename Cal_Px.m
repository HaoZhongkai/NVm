function Px = Cal_Px(t,wh,th)
global Ta Tb N
%%%输入参数为2时意味着wh为AB
switch nargin
    case 2
        n = length(wh)/2;
        ker1 = kernal(0,0);
        kernalArray = repmat(ker1,1,n);
        M0 = ones(length(t),1);
        for i = 1:n
            kernalArray(i).A = wh(2*i-1);
            kernalArray(i).B = wh(2*i);
            M0 = M0.*kernalArray(i).get_M(t);
        end
        Px = 1/2*M0.*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb);
    case 3
        A = wh.*cos(th);
        B = wh.*sin(th);
        n = length(A);
        ker1 = kernal(0,0);
        kernalArray = repmat(ker1,1,n);
        M0 = ones(length(t),1);
        for i = 1:n
            kernalArray(i).A = A(i);
            kernalArray(i).B = B(i);
            M0 = M0.*kernalArray(i).get_M(t);
        end
        Px = 1/2*M0;
end