%% Problem Settings
% Cost Function is at the last of the files
%-------------产生数据
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.01;
t = (tstep:tstep:tmax)';
%%衰减周期,误差,脉冲数,中心信号参数,这里不用考虑衰减
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*20;
%核的参数
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
S0 = Kernal(wh,th,wl,N,t);
S0.get_Px();
S0.AddCentralSignal(N_center,wh_center);
%% 初始化随机背景信号40个核,10组
noise_num = 5;
wh_c = wh_center*rand(40,noise_num);
th_c = pi*rand(40,noise_num);
global M_center;
M_center = zeros(length(t),noise_num);
for i = 1:noise_num
    M_center(:,i) = 2*Get_Px(t,wh_c(:,i),th_c(:,i),wl,N)-1;
end


nVar=12;                % Number of Unknown (Decision) Variables

VarSize=[1 nVar];       % Decision Variables Matrix Size

VarMin=-80e-3;             % Lower Bound of Decision Variables
VarMax= 80e-3;             % Upper Bound of Decision Variables

%% CMA-ES Settings

% Maximum Number of Iterations
MaxIt=300;

% Population Size (and Number of Offsprings)
lambda=(4+round(3*log(nVar)))*10;
% lambda = 300;

% Number of Parents
mu=round(lambda/2);

% Parent Weights
w=log(mu+0.5)-log(1:mu);
w=w/sum(w);

% Number of Effective Solutions
mu_eff=1/sum(w.^2);

% Step Size Control Parameters (c_sigma and d_sigma);
sigma0=0.3*(VarMax-VarMin);
cs=(mu_eff+2)/(nVar+mu_eff+5);
ds=1+cs+2*max(sqrt((mu_eff-1)/(nVar+1))-1,0);
ENN=sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

% Covariance Update Parameters
cc=(4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
c1=2/((nVar+1.3)^2+mu_eff);
alpha_mu=2;
cmu=min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
hth=(1.4+2/(nVar+1))*ENN;

%% Initialization

ps=cell(MaxIt,1);
pc=cell(MaxIt,1);
C=cell(MaxIt,1);
sigma=cell(MaxIt,1);

ps{1}=zeros(VarSize);
pc{1}=zeros(VarSize);
C{1}=eye(nVar);
sigma{1}=sigma0;

empty_individual.Position=[];
empty_individual.Step=[];
empty_individual.Cost=[];

M=repmat(empty_individual,MaxIt,1);
M(1).Position=unifrnd(VarMin,VarMax,VarSize);
M(1).Step=zeros(VarSize);
M(1).Cost=CostFunction(M(1).Position,S0);

BestSol=M(1);

BestCost=zeros(MaxIt,1);

%% CMA-ES Main Loop

for g=1:MaxIt
    
    % Generate Samples
    pop=repmat(empty_individual,lambda,1);
    for i=1:lambda
        pop(i).Step=mvnrnd(zeros(VarSize),C{g});
        pop(i).Position=M(g).Position+sigma{g}*pop(i).Step;
        pop(i).Cost=CostFunction(pop(i).Position,S0);
        
        % Update Best Solution Ever Found
        if pop(i).Cost<BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
  
    % Save Results
    BestCost(g)=BestSol.Cost;
    
    % Display Results
    disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(BestCost(g))]);
    
    % Exit At Last Iteration
    if g==MaxIt
        break;
    end
        
    % Update Mean
    M(g+1).Step=0;
    for j=1:mu
        M(g+1).Step=M(g+1).Step+w(j)*pop(j).Step;
    end
    M(g+1).Position=M(g).Position+sigma{g}*M(g+1).Step;
    M(g+1).Cost=CostFunction(M(g+1).Position,S0);
    if M(g+1).Cost<BestSol.Cost
        BestSol=M(g+1);
    end
    
    % Update Step Size
    ps{g+1}=(1-cs)*ps{g}+sqrt(cs*(2-cs)*mu_eff)*M(g+1).Step/chol(C{g})';
    sigma{g+1}=sigma{g}*exp(cs/ds*(norm(ps{g+1})/ENN-1))^0.3;
    
    % Update Covariance Matrix
    if norm(ps{g+1})/sqrt(1-(1-cs)^(2*(g+1)))<hth
        hs=1;
    else
        hs=0;
    end
    delta=(1-hs)*cc*(2-cc);
    pc{g+1}=(1-cc)*pc{g}+hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
    C{g+1}=(1-c1-cmu)*C{g}+c1*(pc{g+1}'*pc{g+1}+delta*C{g});
    for j=1:mu
        C{g+1}=C{g+1}+cmu*w(j)*pop(j).Step'*pop(j).Step;
    end
    
    % If Covariance Matrix is not Positive Defenite or Near Singular
    [V, E]=eig(C{g+1});
    if any(diag(E)<0)
        E=max(E,0);
        C{g+1}=V*E/V;
    end
    
end

%% Display Results

figure;
% plot(BestCost, 'LineWidth', 2);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

A_best = BestSol.Position(1:end/2);
B_best = BestSol.Position(end/2+1:end);
best_Px = Get_Px2(t,A_best,B_best,wl,N);
figure
plot(t,best_Px,t,S0.Px);
figure
hold on;
scatter(A_best,B_best,'o','SizeData',20);
scatter(A0,B0,'Marker','o','SizeData',30);

%直接使用MSE
function loss = CostFunction(Position,S0)
    global M_center
    A = Position(1:length(Position)/2);
    B = Position(length(Position)/2+1:end);
    noise = M_center(:,randi(5));
    M = 2*Get_Px2(S0.t,A,B,S0.wl,S0.N)-1;
    Px = (M.*noise+1)/2;
    delta = Px-S0.Px;
    loss = delta'*delta;
end
 
