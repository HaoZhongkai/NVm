%ʹ��matlab������
%��һ���汾,������������,����SGD����,sigmoid�����
classdef Network0<handle
    properties
        size        %������ṹ
        num_layers  %���������
        biases      %ƫ��(����ֵ)
        weights     %Ȩ��
    end
    methods
        
        %��ʼ��
        function self = Network0(size0)
            self.size = size0;
            self.num_layers = length(size0);
            %��ʼ��Ȩ����ƫ��(��ֵ)���洢��Ԫ��������
            self.biases = cell(1,length(size0)-1);
            self.weights = cell(1,length(size0)-1);
            for i = 1:(length(size0)-1)
                self.weights{i} = rand(size0(i),size0(i+1));
                self.biases{i} = rand(size0(i),1);
            end
        end
        
        %sigmoid�����
        function z = sigmoid(z0)
            z = 1./(1+exp(-z0));
        end
        
        %sigmoid�����ĵ���
        function z = sigmoid_prime(z0)
            z = sigmoid(z0)*(1-sigmoid(z0));
        end
        
        %ǰ�򴫲�(����Ϊ��һ��Ľ���͸ò�����,��Ϊ������)
        %����ʱע���i���Ȩ�أ�ƫ��洢�ڵ�i-1��Ԫ����
        function output = single_layer_forward(self,input,layer)
            output = input*self.weights{layer-1}+self.biases{layer-1};
            output = sigmoid(output);        
        end
        
    end
end