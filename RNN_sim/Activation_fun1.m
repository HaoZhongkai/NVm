classdef Activation_fun1<handle
    properties
        fun
        param
    end
    methods
        function self = Activation_fun1(at_fun,param)
            switch at_fun
                case 'SR'
                    self.fun = @SymmetricRamp;
                case 'TANH'
                    self.fun = @TANH;
                case 'Sigmoid'
                    self.fun = @Sigmoid;
            
            end
            self.param = param;
        end
        
        function output = activate(self,input)
            output = self.fun(self,input);
        end
            %分别为上界与下界
        function output = SymmetricRamp(self,input)
            inf0 = self.param(1);
            sup0 = self.param(2);
            output = arrayfun(@(x) min(max(x,inf0),sup0),input);
        end
            %a控制tanh的收缩与扩展,b控制幅度
        function output = TANH(self,input)
            a = self.param(1);
            b = self.param(2);
            output = b*tanh(a*input);
        end
        %a控制x方向,b调节纵轴输出
        function output = Sigmoid(self,input)
            a = self.param(1);
            b = self.param(2);
            output = b./(1+exp(-a*input));
        end
    end
    
end