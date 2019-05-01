function J = Jacobi2(t,AB)
step = 1e-5;
J = zeros(length(t),length(AB));
for i = 1:length(AB)
    AB1 = AB;
    AB1(i) = AB(i)+step;
    J(:,i) = (Cal_Px(t,AB1)-Cal_Px(t,AB))/step;
end