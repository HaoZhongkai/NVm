function J = JacobiMatrix(t,AB)
global Ta N
J = zeros(length(t),length(AB));
Px_der = zeros(length(t),length(AB));
n = length(AB)/2;
ker1 = kernal(0,0);
kernalArray = repmat(ker1,1,n);
M = ones(length(t),n);
M0 = ones(length(t),1);
for i = 1:n
    kernalArray(i).A = AB(2*i-1);
    kernalArray(i).B = AB(2*i);
    Px_der(:,2*i-1:2*i) = 1/2*get_M_der(kernalArray(i),t).*exp(-2*N*t/Ta);
    M(:,i) = kernalArray(i).get_M(t);
    M0 = M0.*M(:,i);
end
for i = 1:n
    J(:,2*i-1) = M0./M(:,i).*Px_der(:,2*i-1);
    J(:,2*i) = M0./M(:,i).*Px_der(:,2*i);
end
AB0 = [0.0782,0.0300,0.0407,0.0235,0.0323,0.0445,-0.0130,0.0139,-0.0221...
    0.0245,0.0158,0.0195];
AB0*2;