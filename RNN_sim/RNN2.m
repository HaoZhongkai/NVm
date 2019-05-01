%�ڶ��汾RNN,�������Լ���������Ż�����
classdef RNN2<handle
    properties
        neuron_num
        biase             %���������洢
        weight            %�þ���洢
        at_fun      %��������
    end
    methods
        %�����ʼ��
        function self = RNN2(param_num,activate_fun)
            %b_maxΪ��ʼ������ʱ��Ȩ�����ֵ
            self.biase = zeros(param_num,1);
            %��ʼ��Ȩ��
            self.weight = zeros(param_num);
            %�����
            self.at_fun = activate_fun;
        end
        
        %input ����Ϊ�����������ͬ��������
        function output = feedforward(self,input,tstep)
            %����ĸ��¹�ʽ
%             output = input + tstep*self.biase + self.weight*input
            output = input+tstep*self.biase;
        end
        
        function Renew(self,biase)
            self.biase = biase;
        end
        
        function output = Activate(self,input)
            output = self.at_fun.activate(input);
        end
            
        
    end
end