%使用matlab神经网络
%第一个版本,单隐层神经网络,采用SGD方法,sigmoid激活函数
classdef Network0<handle
    properties
        size        %神经网络结构
        num_layers  %神经网络层数
        biases      %偏差(即阈值)
        weights     %权重
    end
    methods
        
        %初始化
        function self = Network0(size0)
            self.size = size0;
            self.num_layers = length(size0);
            %初始化权重与偏差(阈值)，存储在元胞数组中
            self.biases = cell(1,length(size0)-1);
            self.weights = cell(1,length(size0)-1);
            for i = 1:(length(size0)-1)
                self.weights{i} = rand(size0(i),size0(i+1));
                self.biases{i} = rand(size0(i),1);
            end
        end
        
        %sigmoid激活函数
        function z = sigmoid(z0)
            z = 1./(1+exp(-z0));
        end
        
        %sigmoid函数的导数
        function z = sigmoid_prime(z0)
            z = sigmoid(z0)*(1-sigmoid(z0));
        end
        
        %前向传播(输入为上一层的结果和该层的序号,均为行向量)
        %计算时注意第i层的权重，偏差存储在第i-1个元素中
        function output = single_layer_forward(self,input,layer)
            output = input*self.weights{layer-1}+self.biases{layer-1};
            output = sigmoid(output);        
        end
        
    end
end