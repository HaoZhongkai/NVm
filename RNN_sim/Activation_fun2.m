classdef Activation_fun2<handle
    properties
        fun
        param
    end
    methods
        function self = Activation_fun2(at_fun,param)
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
            %�ֱ�Ϊ�Ͻ����½�
        function output = SymmetricRamp(self,input)
            len = length(input)/2;
            num0 = self.param(1);
            sup0 = self.param(2);
            sup1 = self.param(3);
            output = zeros(size(input));
            index1 = [1:num0,len+1:len+num0];
            index2 = [num0+1:len,len+num0+1:2*len];
            output(index1) = arrayfun(@(x) min(max(x,-sup0)...
                ,sup0),input(index1));
            output(index2) = arrayfun(@(x) min(max(x,...
                -sup1),sup1),input(index2));
        end
            %a����tanh����������չ,b���Ʒ���
        function output = TANH(self,input)
            a = self.param(1);
            b = self.param(2);
            output = b*tanh(a*input);
        end
        %a����x����,b�����������
        function output = Sigmoid(self,input)
            a = self.param(1);
            b = self.param(2);
            output = b./(1+exp(-a*input));
        end
    end
    
end