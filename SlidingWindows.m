%-----------------
%这个函数为分离尖峰的算法，输入应为信号函数的导数
%分离出的baseinfo为基线的下标列表
%peakinfo为所有峰的结构数组,包含起止点
%若Nopeak为1，则没有找到峰
%记得Px是一个列向量
function [BaseIndex,PeakInfo,NoPeak] = SlidingWindows(DiffData,NoiseFactor)
    DiffDataSize = length(DiffData);
    DiffData = abs(DiffData);
    HalfWindow = 0.015*DiffDataSize;



    %---------
    %先求出误差的标准差
    NoisePoint = DiffDataSize / 20;
    tempNoise = zeros(1,20);
    for i = 1:20
        tempNoise(i) = std(DiffData((i-1)*NoisePoint+1:i*NoisePoint));
    end
    Noise = min(tempNoise);
    Noise = Noise * NoiseFactor;


    %--------------
    %类似于神经网络的padding技巧，两侧补数据
    padding_data = zeros(HalfWindow,1);
    DiffData = [padding_data;DiffData;padding_data];
    [TempWindowMax,TempMaxIndex] = max(DiffData(1:2*HalfWindow+1));
    [TempWindowMin,TempMinIndex] = min(DiffData(1:2*HalfWindow+1));
    %用于指示上一个点是否是峰的状态变量
    PeakState = 0;
    PeakStart = 0;
    PeakNum = 0;
    %-----------
    %开始寻峰
    %未了避免对每个点都重新算一次窗高，采用局部调整的方式来计算
    for i = 1:DiffDataSize
        %调整当前窗的最大值
        if(TempMaxIndex>=i-HalfWindow && DiffData(i+HalfWindow)>TempWindowMax)
            TempWindowMax = DiffData(i+HalfWindow);
            TempMaxIndex = i+HalfWindow;
        elseif(i-HalfWindow>TempMaxIndex)
            [TempWindowMax,TempMaxIndex] = max(DiffData(i-HalfWindow:i+HalfWindow)); 
        end

        %调整当前窗的最小值
        if(TempMinIndex>=i-HalfWindow && DiffData(i+HalfWindow)<TempWindowMin)
            TempWindowMin = DiffData(i+HalfWindow);
            TempMinIndex = i+HalfWindow;
        elseif(i-HalfWindow>TempMinIndex)
            [TempWindowMin,TempMinIndex] = min(DiffData(i-HalfWindow:i+HalfWindow)); 
        end
        
        %判断是否为峰
        TempWindowHeight = TempWindowMax-TempWindowMin;
        if(TempWindowHeight>=Noise)
            %该点为峰但前一个点不为峰的时候,记下初始点
            if(PeakState==0)
                PeakStart = i;
                PeakState = 1;
            end
        elseif(TempWindowHeight<=Noise)
            %该点不为峰但前一个点为峰时，记下峰的信息
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


