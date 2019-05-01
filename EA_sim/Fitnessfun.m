%损失函数类
classdef Fitnessfun<handle
    properties
        flag
        param
        fun
    end
    methods
        function self = Fitnessfun(flag,param)
            self.flag = flag;
            self.param = param;
            switch flag
                case 'MSE'
                    self.fun = @MSEfitness;
                case 'Correlation'
                    self.fun = @Correlation;
                case 'MIX'
                    self.fun  = @MIXfitness;
            end
        end
        
        function fitness = Get_fit(self,S1,S0)
            fitness = self.fun(self,S1,S0);
        end
        
        function fitness = MSEfitness(self,S1,S0)
            %MSEfitness是MSEloss的倒数乘一个值
            ratio = self.param;
            P1 = 1-S1.Px;
            P0 = 1-S0.Px;
            P1(P1<0.05)=0;
            P0(P0<0.05)=0;
            fitness = ratio/((P1-P0)'*(P1-P0)*(S0.t(2)...
                -S0.t(1)));
        end
        %关联度,在0处乘起来积分
        function fitness = Correlation(self,S1,S0)
            ratio = self.param;
            P1 = S1.Px;
            P0 = S0.Px;
            fitness = P1'*P0/(sqrt(P1'*P1)*ratio);
        end
        %将交叉关联函数归一化再除以归一化的sqrt(MSEloss)
        function fitness = MIXfitness(self,S1,S0)
            ratio = self.param;
            P1 = 1-S1.Px;
            P0 = 1-S0.Px; 
            P1(P1<0.1)=0;
            P0(P0<0.1)=0;
%             fitness = P1'*P0/sqrt((P1'*P1)*exp(ratio*(P1-P0))'*...
%                 exp(ratio*(P1-P0)));
            fitness = P1'*P0/sqrt((P1'*P1)*(P1-P0)'*(P1-P0));
        end
    end
end