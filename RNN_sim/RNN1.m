%RNN第一版本,仅含单层神经元
%因为参数个数都是偶数个,所以输入的参数个数为节点个数的1/2
classdef RNN1<handle
    properties
        neuron_num
        biase             %用列向量存储
        weight            %用矩阵存储
        at_fun      %激活函数句柄
    end
    methods
        %网络初始化
        function self = RNN1(param_num,activate_fun)
            %b_max为初始化网络时的权重最大值
            self.biase = zeros(2*param_num,1);
            %初始化权重
            self.weight = zeros(param_num);
            %激活函数
            self.at_fun = activate_fun;
        end
        
        %input 必须为与参数个数相同的列向量
        function output = feedforward(self,input,tstep)
            %暂时因为weight==0,暂时不考虑
            output = input - tstep*self.biase;
        end
        
        function Renew(self,biase)
            self.biase = biase;
        end
        
        function output = Activate(self,input)
            output = self.at_fun.activate(input);
        end
            
        
    end
end