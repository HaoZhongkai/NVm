%-------------
%本程序为将峰(未做对多核的处理)进行等距分离的示例程序，可以直接运行
%该函数在单一核的峰划分上取得了比较好的效果
%考虑到多核时数据结构会变复杂，将直接用python重构本函数
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
tmax = 10;
tstep = 0.01;
t = (tstep:tstep:tmax)';
Ta = 5e3;
Tb = 5e3;
%set parameter...
wh = 1e-3*[83.8];
th = pi/180*[21];
N = 56; 
wl = 2*pi*gama0*B0;
wb = sqrt((2*pi*wh).^2+2*2*pi*wl*wh.*cos(th)+wl^2);
A = wh.*cos(th);
B = wh.*sin(th);
Px = Get_Px(t,wh,th,wl,N);
%第一个中心峰位置
T0 = pi/(2*wl);
%任意给定一个点的下标,然后它将自动在周期中搜索
%搜索所有与中心峰周期数相同的峰
N_step = 5;
N_near = 25;
%搜索中心峰周围N00个点，步长为N_step
% N_center = floor(T0/tstep);
N0 = 532;
period_num = floor((length(Px)/N0-1)/2)+1;
N00 = 0;
N_list = N0-N00:N_step:N0+N00;
signal = zeros(length(N_list),period_num);
index_signal = zeros(length(N_list),period_num);
t_signal = zeros(length(N_list),period_num);
for i = 1:length(N_list)
    [index0, signal0] = Period_search(t,Px,N_list(i),N_near);
    index_signal(i,:) = index0;
    t_signal(i,:) = t(index0)';
    signal(i,:) = signal0;
end
% figure
% plot(t,Px,'color','b');
Px2 = 1-Px;
%-----------------------
%使用matlab寻峰算法寻峰
[height,loc,pk_w,pk_p] = findpeaks(Px2,t,'MinPeakProminence',0.03,...
    'MinPeakDistance',0.05,'Annotate','extents');
% for i = 1:length(N_list)
%     s = scatter(t_signal(i,:),signal(i,:),'MarkerFaceColor',rand(1,3));
% end
figure
scatter(loc,height,'MarkerEdgeColor','r');
% plot(t,Px2,'color','b');

%-------------------
%求出差矩阵
loc_diff = diff_peaklocs(loc);
% for i = 1:length(loc)
%     scatter(i+1:length(loc),loc_diff(i,i+1:end));
% end


%-----------------
%从中挑出一个间距，并类推然后求出匹配程度

%以第一行为例
%将loc换成行向量
match_threshold = 1e-5;     %匹配时跳出的阈值
loc = loc';
pk_len = length(loc);
loc_max = loc(end);
S = zeros(1,floor(1/3*pk_len));
R2 = S;
%注意S已换算成比值
for i = 1:floor(1/3*pk_len)
    pk_d = loc_diff(1,i+1);
    npk = floor((loc_max-loc(1))/pk_d)+1;
    loc2 = loc(1)+(0:npk-1)*pk_d;
    [loc_m,loc_index] = pk_match(loc2,loc(1,:));
    [loc2,loc_m,~] = outlier_remove(loc2,loc_m,loc_index,3.5);
    S(i) = L2biase(loc2,loc_m);
    R2(i) = Rsquare(loc2,loc_m);
    [pk_h,pk_index] = matchpeak(Px,loc_m,tstep);
    if R2(i)>1-match_threshold
        break;
    end
%     画所有峰与Px上的图
%     figure
%     hold on; 
%     plot(t,Px);
%     scatter(t(pk_index),pk_h);
%     textstr = ['i = ',num2str(i),' S = ',num2str(S(i)),' R2 = ',num2str(R2(i))];
%     annotation('textbox',[.2 .5 .3 .3],'String',textstr,'FitBoxtoText','on');

%     画匹配出的峰的位置序列在所有峰的位置序列中的位置
%     figure
%     hold on;
%     scatter(1:length(loc(1,:)),loc(1,:));
%     scatter(loc_index,loc_m,'+');

%     画理论的峰的序列与数据匹配程度的图
%     figure
%     hold on;
%     scatter(1:length(loc_m),loc_m,'o');
%     plot(1:length(loc_m),loc2,'color','b');
end
[Smin,best_match_index] = min(S);
[loc_m,loc_index] = pk_match(loc2,loc(1,:));
[loc2,loc_m,loc_index] = outlier_remove(loc2,loc_m,loc_index,3.5);
[pk_h,pk_index] = matchpeak(Px,loc_m,tstep);


%     画理论的峰的序列与数据匹配程度的图
figure
hold on;
scatter(1:length(loc_index),loc(loc_index));
scatter(1:length(loc_index),loc2,'+');

%画所有峰与Px上的图
figure
hold on;
plot(t,Px);
scatter(t(pk_index),pk_h);
textstr = ['i = ',num2str(best_match_index),' S = ',num2str(S(...
    best_match_index)),' R2 = ',num2str(R2(best_match_index))];
annotation('textbox',[.2 .5 .3 .3],'String',textstr,'FitBoxtoText','on');

%对峰的位置中的任意两个点都求差，这样得出一个矩阵
%y(i,j) = y(j) - y(i)
function loc_diff = diff_peaklocs(pk_locs)
    loc_diff = zeros(length(pk_locs));
    for i = 1:length(pk_locs)
        loc_diff(i,:) = pk_locs(:)-pk_locs(i);
    end
end



%----------------
%找到y0最匹配的序列
%最后匹配到的序列可以重复(考虑到可能会有峰缺失的情况)
function [y_match,index] = pk_match(y,y0)
    pivot = 1;
    y_m = 0;
    y_match = zeros(size(y));
    index = zeros(1,length(y));
    for i = 1:length(y)
        %寻找匹配点
        while pivot<=length(y0)&&abs(y_m-y(i))>=abs(y(i)-y0(pivot))
            y_m = y0(pivot);
            pivot = pivot+1;
        end
        y_match(i) = y_m;
        index(i) = pivot-1;
        pivot = max(1,pivot-1);       %确保所有点都能匹配上(暂不考虑准确性)
    end
end


%输入y为算出的个峰的位置,y0为数据(已经经过匹配),计算偏差
%暂时以L2范数为例
function S = L2biase(y,y0)
%     S = sqrt((y-y0)*(y-y0)')/(length(y)*(y(2)-y(1)));
    S = sum(abs(y-y0))/(length(y)*(y(2)-y(1)));
end

%计算相关系数,下标为1到n
%y0为数据点，y为拟合直线,
function R2 = Rsquare(y,y0)
    y0_b = sum(y0)/length(y0);
    R2 = 1-((y0-y)*(y0-y)')/((y0-y0_b)*(y0-y0_b)');
end


%--------------------
%使用modified Z-score方法来剔除异常点
%该序列已经经过匹配，为等长序列
%y0_index为输出的向量在所有峰数据中的下标索引
function [y_out,y0_out,y0_index] = outlier_remove(y,y0,index,threshold)
    MAD = median(abs(y0-y));
    d_median = median(y0-y);
    param = 0.6745;
    i = 1;
    while i<=length(y)
        if param*abs(y0(i)-y(i)-d_median)/MAD>threshold
            y0(i) = [];
            y(i) = [];
            index(i) = [];
        end
        i = i+1;
    end
    y0_index = index;
    y_out = y;
    y0_out = y0;
end

%------------------
%在Px中匹配该峰的位置
%tstep为t的步长，为了匹配t
function [pk_h,pk_index] = matchpeak(Px,loc,tstep)
    pk_index = round(loc/tstep);
    pk_h = Px(pk_index);
end
    
    
    

