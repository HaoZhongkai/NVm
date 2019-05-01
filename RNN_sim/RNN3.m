%第二版本RNN,解决加上约束条件的优化问题
classdef RNN3<handle
    properties
        neuron_num
        biase             %用列向量存储
        weight            %用矩阵存储
        at_fun      %激活函数句柄
    end
    methods
        %网络初始化
        function self = RNN3(param_num,activate_fun)
            %b_max为初始化网络时的权重最大值
            self.biase = zeros(param_num,1);
            %初始化权重
            self.weight = zeros(param_num);
            %激活函数
            self.at_fun = activate_fun;
        end
        
        %input 必须为与参数个数相同的列向量
        function output = feedforward(self,input,tstep)
            %输出的更新公式
            output = input + tstep*self.biase + self.weight*input;
%             output = input+tstep*self.biase;
        end
        
        function Renew(self,weight,biase)
            self.biase = biase;
            self.weight = weight;
        end
        
        function output = Activate(self,input)
            output = self.at_fun.activate(input);
        end
            
        
    end
end