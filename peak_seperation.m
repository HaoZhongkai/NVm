%-------------
%������Ϊ����(δ���Զ�˵Ĵ���)���еȾ�����ʾ�����򣬿���ֱ������
%�ú����ڵ�һ�˵ķ廮����ȡ���˱ȽϺõ�Ч��
%���ǵ����ʱ���ݽṹ��临�ӣ���ֱ����python�ع�������
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
%��һ�����ķ�λ��
T0 = pi/(2*wl);
%�������һ������±�,Ȼ�������Զ�������������
%�������������ķ���������ͬ�ķ�
N_step = 5;
N_near = 25;
%�������ķ���ΧN00���㣬����ΪN_step
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
%ʹ��matlabѰ���㷨Ѱ��
[height,loc,pk_w,pk_p] = findpeaks(Px2,t,'MinPeakProminence',0.03,...
    'MinPeakDistance',0.05,'Annotate','extents');
% for i = 1:length(N_list)
%     s = scatter(t_signal(i,:),signal(i,:),'MarkerFaceColor',rand(1,3));
% end
figure
scatter(loc,height,'MarkerEdgeColor','r');
% plot(t,Px2,'color','b');

%-------------------
%��������
loc_diff = diff_peaklocs(loc);
% for i = 1:length(loc)
%     scatter(i+1:length(loc),loc_diff(i,i+1:end));
% end


%-----------------
%��������һ����࣬������Ȼ�����ƥ��̶�

%�Ե�һ��Ϊ��
%��loc����������
match_threshold = 1e-5;     %ƥ��ʱ��������ֵ
loc = loc';
pk_len = length(loc);
loc_max = loc(end);
S = zeros(1,floor(1/3*pk_len));
R2 = S;
%ע��S�ѻ���ɱ�ֵ
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
%     �����з���Px�ϵ�ͼ
%     figure
%     hold on; 
%     plot(t,Px);
%     scatter(t(pk_index),pk_h);
%     textstr = ['i = ',num2str(i),' S = ',num2str(S(i)),' R2 = ',num2str(R2(i))];
%     annotation('textbox',[.2 .5 .3 .3],'String',textstr,'FitBoxtoText','on');

%     ��ƥ����ķ��λ�����������з��λ�������е�λ��
%     figure
%     hold on;
%     scatter(1:length(loc(1,:)),loc(1,:));
%     scatter(loc_index,loc_m,'+');

%     �����۵ķ������������ƥ��̶ȵ�ͼ
%     figure
%     hold on;
%     scatter(1:length(loc_m),loc_m,'o');
%     plot(1:length(loc_m),loc2,'color','b');
end
[Smin,best_match_index] = min(S);
[loc_m,loc_index] = pk_match(loc2,loc(1,:));
[loc2,loc_m,loc_index] = outlier_remove(loc2,loc_m,loc_index,3.5);
[pk_h,pk_index] = matchpeak(Px,loc_m,tstep);


%     �����۵ķ������������ƥ��̶ȵ�ͼ
figure
hold on;
scatter(1:length(loc_index),loc(loc_index));
scatter(1:length(loc_index),loc2,'+');

%�����з���Px�ϵ�ͼ
figure
hold on;
plot(t,Px);
scatter(t(pk_index),pk_h);
textstr = ['i = ',num2str(best_match_index),' S = ',num2str(S(...
    best_match_index)),' R2 = ',num2str(R2(best_match_index))];
annotation('textbox',[.2 .5 .3 .3],'String',textstr,'FitBoxtoText','on');

%�Է��λ���е����������㶼�������ó�һ������
%y(i,j) = y(j) - y(i)
function loc_diff = diff_peaklocs(pk_locs)
    loc_diff = zeros(length(pk_locs));
    for i = 1:length(pk_locs)
        loc_diff(i,:) = pk_locs(:)-pk_locs(i);
    end
end



%----------------
%�ҵ�y0��ƥ�������
%���ƥ�䵽�����п����ظ�(���ǵ����ܻ��з�ȱʧ�����)
function [y_match,index] = pk_match(y,y0)
    pivot = 1;
    y_m = 0;
    y_match = zeros(size(y));
    index = zeros(1,length(y));
    for i = 1:length(y)
        %Ѱ��ƥ���
        while pivot<=length(y0)&&abs(y_m-y(i))>=abs(y(i)-y0(pivot))
            y_m = y0(pivot);
            pivot = pivot+1;
        end
        y_match(i) = y_m;
        index(i) = pivot-1;
        pivot = max(1,pivot-1);       %ȷ�����е㶼��ƥ����(�ݲ�����׼ȷ��)
    end
end


%����yΪ����ĸ����λ��,y0Ϊ����(�Ѿ�����ƥ��),����ƫ��
%��ʱ��L2����Ϊ��
function S = L2biase(y,y0)
%     S = sqrt((y-y0)*(y-y0)')/(length(y)*(y(2)-y(1)));
    S = sum(abs(y-y0))/(length(y)*(y(2)-y(1)));
end

%�������ϵ��,�±�Ϊ1��n
%y0Ϊ���ݵ㣬yΪ���ֱ��,
function R2 = Rsquare(y,y0)
    y0_b = sum(y0)/length(y0);
    R2 = 1-((y0-y)*(y0-y)')/((y0-y0_b)*(y0-y0_b)');
end


%--------------------
%ʹ��modified Z-score�������޳��쳣��
%�������Ѿ�����ƥ�䣬Ϊ�ȳ�����
%y0_indexΪ��������������з������е��±�����
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
%��Px��ƥ��÷��λ��
%tstepΪt�Ĳ�����Ϊ��ƥ��t
function [pk_h,pk_index] = matchpeak(Px,loc,tstep)
    pk_index = round(loc/tstep);
    pk_h = Px(pk_index);
end
    
    
    

