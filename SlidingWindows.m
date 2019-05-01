%-----------------
%�������Ϊ��������㷨������ӦΪ�źź����ĵ���
%�������baseinfoΪ���ߵ��±��б�
%peakinfoΪ���з�Ľṹ����,������ֹ��
%��NopeakΪ1����û���ҵ���
%�ǵ�Px��һ��������
function [BaseIndex,PeakInfo,NoPeak] = SlidingWindows(DiffData,NoiseFactor)
    DiffDataSize = length(DiffData);
    DiffData = abs(DiffData);
    HalfWindow = 0.015*DiffDataSize;



    %---------
    %��������ı�׼��
    NoisePoint = DiffDataSize / 20;
    tempNoise = zeros(1,20);
    for i = 1:20
        tempNoise(i) = std(DiffData((i-1)*NoisePoint+1:i*NoisePoint));
    end
    Noise = min(tempNoise);
    Noise = Noise * NoiseFactor;


    %--------------
    %�������������padding���ɣ����ಹ����
    padding_data = zeros(HalfWindow,1);
    DiffData = [padding_data;DiffData;padding_data];
    [TempWindowMax,TempMaxIndex] = max(DiffData(1:2*HalfWindow+1));
    [TempWindowMin,TempMinIndex] = min(DiffData(1:2*HalfWindow+1));
    %����ָʾ��һ�����Ƿ��Ƿ��״̬����
    PeakState = 0;
    PeakStart = 0;
    PeakNum = 0;
    %-----------
    %��ʼѰ��
    %δ�˱����ÿ���㶼������һ�δ��ߣ����þֲ������ķ�ʽ������
    for i = 1:DiffDataSize
        %������ǰ�������ֵ
        if(TempMaxIndex>=i-HalfWindow && DiffData(i+HalfWindow)>TempWindowMax)
            TempWindowMax = DiffData(i+HalfWindow);
            TempMaxIndex = i+HalfWindow;
        elseif(i-HalfWindow>TempMaxIndex)
            [TempWindowMax,TempMaxIndex] = max(DiffData(i-HalfWindow:i+HalfWindow)); 
        end

        %������ǰ������Сֵ
        if(TempMinIndex>=i-HalfWindow && DiffData(i+HalfWindow)<TempWindowMin)
            TempWindowMin = DiffData(i+HalfWindow);
            TempMinIndex = i+HalfWindow;
        elseif(i-HalfWindow>TempMinIndex)
            [TempWindowMin,TempMinIndex] = min(DiffData(i-HalfWindow:i+HalfWindow)); 
        end
        
        %�ж��Ƿ�Ϊ��
        TempWindowHeight = TempWindowMax-TempWindowMin;
        if(TempWindowHeight>=Noise)
            %�õ�Ϊ�嵫ǰһ���㲻Ϊ���ʱ��,���³�ʼ��
            if(PeakState==0)
                PeakStart = i;
                PeakState = 1;
            end
        elseif(TempWindowHeight<=Noise)
            %�õ㲻Ϊ�嵫ǰһ����Ϊ��ʱ�����·����Ϣ
            if(PeakState==1)
                PeakEnd = i-1;
                PeakState = 0;
                PeakNum = PeakNum+1;
                PeakInfo(PeakNum) = struct('Start',PeakStart,'End',PeakEnd);
            end
        end
    end        


    BaseIndex = 1:DiffDataSize;
    NoPeak = 0;
    if(PeakNum == 0)
        disp('you have find no peak in this spectra');
        PeakInfo = struct([]);
        NoPeak=1;
    else
        DelIndex = cell(1,length(PeakInfo));
        for j = 1:length(PeakInfo)
            DelIndex{j} = PeakInfo(j).Start:PeakInfo(j).End; 
        end
        BaseIndex(horzcat(DelIndex{:})) = [];
    end
end


