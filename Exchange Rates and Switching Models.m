%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modification of code for ex 5 

clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Import Data from file     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sheetName='Foglio1';
[numbers, strings] = xlsread('Exchange.xlsx', sheetName);
if ~isempty(numbers)
    newData1.data =  numbers;
end
if ~isempty(strings)
    newData1.textdata =  strings;
end
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Transform Data            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sm=0;    % window for smoothing data, 0 do not smooth data
Y=data(:,1); % 1: 'EUUS' 2: 'YEUS' 3:'UKUS'
YY=(log(Y(2:end))-log(Y(1:end-1)))*100;
T=size(YY,1);
Z=data(:,2); % 1: 'EUUS' 2: 'YEUS' 3:'UKUS'
ZZ=(log(Z(2:end))-log(Z(1:end-1)))*100;

if sm>0
    X(1:sm-1,:)=YY(1:sm-1,1);
    for i=sm:T
        X(i,:)=mean(YY(i-sm+1:i,1));
    end
else
    X=YY;
end
if sm>0
    Z(1:sm-1,:)=ZZ(1:sm-1,1);
    for j=sm:T
        Z(j,:)=mean(ZZ(j-sm+1:j,1));
    end
else
    Z_1=ZZ;
end
startDate = datenum('01-01-1999');
endDate=datenum('11-03-2011');
xData = linspace(startDate,endDate,T);
zLAG=lagmatrix(Z_1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Inference - Switching Model   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
M=2000;                  %total number of MCMC iterations
burn=floor(M*0.5);   %proportion of initial burn-in MCMC iterations

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Prior setting
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Location prior distribution
m1=0.05;
m0=-0.05;
gam1=0.05;
gam0=0.05;
sigma1=0.05;
sigma0=0.05;
mu1=0.05;
mu0=-0.05;
%%%%%%%%%% Scale prior distribution
a1=10;
b1=4;
a0=10;
b0=4;
%%%%%%%%%% Mixing prob. prior distribution
c=10;
d=30;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
alf1_hat=zeros(M,1);
alf0_hat=zeros(M,1);
theta_hat=zeros(M,1);
sig21_hat=zeros(M,1);
sig20_hat=zeros(M,1);
beta1_hat=zeros(M,1);
beta0_hat=zeros(M,1);
U_hat=zeros(M,T);
sig21_hat(1,1)=0.02;
sig20_hat(1,1)=0.02;
alf1_hat(1,1)=0.1;
alf0_hat(1,1)=-0.1;
beta1_hat(1,1)=0.1;
beta0_hat(1,1)=-0.1;
theta_hat(1,1)=0.1;
%U_hat(1,:)=(rand(T,1)<0.08);
i=1;
for t=1:T
    wmean1_hat=alf1_hat(i,1)-beta1_hat(i,1)*zLAG(t,1);
    wmean0_hat=alf0_hat(i,1)-beta0_hat(i,1)*zLAG(t,1);
    w=theta_hat(i,1)*normpdf(X(t,1),alf1_hat(i,1),sqrt(sig21_hat(i,1)))...
        /(theta_hat(i,1)*normpdf(X(t,1),alf1_hat(i,1),sqrt(sig21_hat(i,1)))+...
        (1-theta_hat(i,1))*normpdf(X(t,1),alf0_hat(i,1),sqrt(sig20_hat(i,1))));
    U_hat(i,t)=rand(1,1)<w;
end
%%
for i=2:M
    %%%%%%%%%%%% full conditional alfa0 and alfa1
    T1=sum(U_hat(i-1,:)==1);
    T0=T-T1;
   
    gam1_bar=(T1/sig21_hat(i-1,1)+1/gam1)^(-1);
    gam0_bar=(T0/sig20_hat(i-1,1)+1/gam0)^(-1);
    m1_bar=gam1_bar*(sum(X((U_hat(i-1,:)==1),1)))-beta1_hat(i,1)*(sum(zLAG((U_hat(i-1,:)==1),1)))/sig21_hat(i-1,1)+m1/gam1;
    m0_bar=gam0_bar*(sum(X(T0,1)))-beta0_hat(i,1)*(sum(zLAG(T0,1)))/sig20_hat(i-1,1)+m0/gam0;
    alf1_hat(i,1)=m1_bar+sqrt(gam1_bar)*random('normal',0,1);
    alf0_hat(i,1)=m0_bar+sqrt(gam0_bar)*random('normal',0,1);
%%%%%%%%%%%%%%%%full conditional on beta0 and beta1 
    sigma0_bar=(1/sig20_hat(i-1,1)*(sum(zLAG(T0,1))^2)+1/sigma0)^(-1);
    sigma1_bar=(1/sig21_hat(i-1,1)*(sum(zLAG((U_hat(i-1,:)==1),1)))^2)+1/sigma1^(-1);
    u0_bar= sigma1_bar*1/sig21_hat(i-1,1)*(sum(X(T0,1))-alf0_hat(i,1)*(sum(zLAG(T0,1))));
    u1_bar= sigma1_bar*(1/sig21_hat(i-1,1)*(sum(X((U_hat(i-1,:)==1),1)))-alf0_hat(i,1)*(sum(zLAG((U_hat(i-1,:)==1),1))));
    beta0_hat(i,1)=u0_bar+sqrt(sigma0_bar)*random('normal',0,1);
    beta1_hat(i,1)=u1_bar+sqrt(sigma1_bar)*random('normal',0,1);
    %%%%%%%%%%%% full conditional sigma2_0 and sigma2_1
    a1_bar=a1+T1/2;
    a0_bar=a0+T0/2;
    e1=((sum(X((U_hat(i-1,:)==1),1)))-alf1_hat(i,1)-(beta1_hat(i,1)*(sum(zLAG((U_hat(i-1,:)==1),1)))));
    e0=(X(T0,1)-alf0_hat(i,1)-(beta0_hat(i,1)*(sum(zLAG(T0,1)))));
    b1_bar=b1+(e1'*e1)/2;
    b0_bar=b0+(e0'*e0)/2;
    sig21_hat(i,1)=1/random('gamma',a1_bar,1/b1_bar);
    sig20_hat(i,1)=1/random('gamma',a0_bar,1/b0_bar);
    while (sig20_hat(i,1)>sig21_hat(i,1))
        sig21_hat(i,1)=1/random('gamma',a1_bar,1/b1_bar);
        sig20_hat(i,1)=1/random('gamma',a0_bar,1/b0_bar);
    end
    %%%%%%%%%%%% full conditional theta
    c_bar=c+T1+1;
    d_bar=d+T0+1;
    theta_hat(i,1)=random('beta',c_bar,d_bar);
    %%%%%%%%%%%% full conditional Ut
    for t=1:T
        w=theta_hat(i,1)*normpdf(X(t,1),alf1_hat(i,1),sqrt(sig21_hat(i,1)))...
            /(theta_hat(i,1)*normpdf(X(t,1),alf1_hat(i,1),sqrt(sig21_hat(i,1)))+...
            (1-theta_hat(i,1))*normpdf(X(t,1),alf0_hat(i,1),sqrt(sig20_hat(i,1))));
        U_hat(i,t)=rand(1,1)<w;
    end
end
% Rolling window estimates
tau=60; % window size
tau1=180;
svol=zeros(T-tau1,2);
for i=tau1+1:T
    svol(i,1)=std(X(i-tau:i,1));
end
for i=tau1+1:T
end 
    svol(i,2)=std(X(i-tau1:i,1));
%%  %%%% Posterior Progressive MCMC averages
figure(13)
subplot(3,2,1)
plot(cumsum(alf1_hat)./(1:M)');%ylim([-0.5 0.5]);
ylabel('$\hat \alpha_{1}$','interpreter','latex')
set(gca,'fontsize',12);
subplot(3,2,2)
plot(cumsum(alf0_hat)./(1:M)');%ylim([-0.5 0.5]);
ylabel('$\hat \alpha_{0}$','interpreter','latex')
set(gca,'fontsize',12);
subplot(3,2,3)
plot(cumsum(beta1_hat)./(1:M)');%ylim([-0.5 0.5]);
ylabel('$\hat \beta1_{0}$','interpreter','latex')
set(gca,'fontsize',12);
plot(cumsum(beta0_hat)./(1:M)');%ylim([-0.5 0.5]);
ylabel('$\hat \beta0_{0}$','interpreter','latex')
set(gca,'fontsize',12);
subplot(3,2,3)
plot(cumsum(sig21_hat)./(1:M)')
ylabel('$\hat \sigma^{2}_{1}$','interpreter','latex')
set(gca,'fontsize',12);
subplot(3,2,4)
plot(cumsum(sig20_hat)./(1:M)')
ylabel('$\hat \sigma^{2}_{0}$','interpreter','latex')
set(gca,'fontsize',12);
subplot(3,2,5:6)
plot(cumsum(theta_hat)./(1:M)')
ylabel('$\hat \theta$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'MCMCav.eps','-depsc');

%%  %%%% Posterior Progressive MCMC trace plots
figure(14)
subplot(3,2,1)
plot(alf1_hat);%ylim([-0.5 0.5]);
xlabel('$i$','interpreter','latex');
ylabel('$\alpha^{(i)}_{1}$','interpreter','latex')
set(gca,'fontsize',12);
%
subplot(3,2,2)
plot(alf0_hat);%ylim([-0.5 0.5]);
xlabel('$i$','interpreter','latex');
ylabel('$\alpha^{(i)}_{0}$','interpreter','latex')
set(gca,'fontsize',12);
%
subplot(3,2,2)
plot(beta1_hat);%ylim([-0.5 0.5]);
xlabel('$i$','interpreter','latex');
ylabel('$\beta1^{(i)}_{0}$','interpreter','latex')
set(gca,'fontsize',12);
%
subplot(3,2,2)
plot(beta0_hat);%ylim([-0.5 0.5]);
xlabel('$i$','interpreter','latex');
ylabel('$\beta0^{(i)}_{0}$','interpreter','latex')
set(gca,'fontsize',12);
%
subplot(3,2,3)
plot(sig21_hat)
xlabel('$i$','interpreter','latex');
ylabel('$\sigma^{2,(i)}_{1}$','interpreter','latex')
set(gca,'fontsize',12);
%
subplot(3,2,4)
plot(sig20_hat)
xlabel('$i$','interpreter','latex');
ylabel('$\sigma^{2,(i)}_{0}$','interpreter','latex')
set(gca,'fontsize',12);
%
subplot(3,2,5:6)
plot(theta_hat)
xlabel('$i$','interpreter','latex');
ylabel('$\theta^{(i)}$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'MCMC.eps','-depsc');

%%  %%%% Posterior MCMC histogram
figure(151)
grd=(-0.1:0.01:0.2);
[f,grd]=histcounts(alf1_hat(burn+1:M,:),grd);
pdfPrior=pdf('normal',grd,m1,sqrt(gam1));%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\alpha_{1}$','interpreter','latex')
ylabel('$p(\alpha_{1})$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'PosteriorAlf1.eps','-depsc');
%
figure(152)
grd=(-0.1:0.01:0.2);
[f,grd]=histcounts(alf0_hat(burn+1:M,:),grd);
pdfPrior=pdf('normal',grd,m0,sqrt(gam0));%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\alpha_{0}$','interpreter','latex')
ylabel('$p(\alpha_{0})$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'PosteriorAlf0.eps','-depsc');
%
figure(152)
grd=(-0.1:0.01:0.2);
[f,grd]=histcounts(beta1_hat(burn+1:M,:),grd);
pdfPrior=pdf('normal',grd,m0,sqrt(gam0));%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\beta_{1}$','interpreter','latex')
ylabel('$p(\beta_{1})$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'PosteriorAlf0.eps','-depsc');
%
figure(152)
grd=(-0.1:0.01:0.2);
[f,grd]=histcounts(beta0_hat(burn+1:M,:),grd);
pdfPrior=pdf('normal',grd,m0,sqrt(gam0));%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\beta_{0}$','interpreter','latex')
ylabel('$p(\beta_{0})$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'PosteriorAlf0.eps','-depsc');
%
figure(153)
grd=(0:0.01:1);
[f,grd]=histcounts(sig21_hat(burn+1:M,:),grd);
pdfPrior=pdfIGamma(grd,a1,b1);%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\sigma^{2}_{1}$','interpreter','latex')
ylabel('$p(\sigma^{2}_{1})$','interpreter','latex')
print(gcf,'PosteriorSig21.eps','-depsc');
%
figure(154)
grd=(0:0.01:1);
[f,grd]=histcounts(sig20_hat(burn+1:M,:),grd);
pdfPrior=pdfIGamma(grd,a0,b0);%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\sigma^{2}_{0}$','interpreter','latex')
ylabel('$p(\sigma^{2}_{0})$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'PosteriorSig20.eps','-depsc');
%
figure(155)
grd=(0:1/100:1);
[f,grd]=histcounts(theta_hat(burn+1:M,:),grd);
pdfPrior=pdf('beta',grd,c,d);%prior pdf
bar(grd(2:end),f/(M-burn));
hold on;
plot(grd,pdfPrior*(grd(2)-grd(1)));
hold off;
xlabel('$\theta$','interpreter','latex')
ylabel('$p(\theta)$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'PosteriorTheta.eps','-depsc');
%%
figure(121)
plot(xData(tau1:end),abs(X(tau1:end)),'color',[0.8 0.8 0.9]);
ylim([0,2]);
hold on;
plot(xData(tau1:end),svol(tau1:end,1),'b--','linewidth',1.2);
plot(xData(tau1:end),svol(tau1:end,2),'b:','linewidth',1.2);
hold off;
ylim([-0.01,2]);
ax = gca;
ax.XTick = xData;
xticks((xData(tau1):800:xData(end)));
datetick('x','dd/mm/yy','keepticks')
legend('$|y_t|$','$\hat{\sigma}_t^{(60)}$','$\hat{\sigma}_t^{(180)}$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'VolHat1.eps','-depsc');

figure(122)
area(xData(tau1:end),(mean(U_hat(burn+1:M,tau1:end))),'facecolor',[0.85,0.8,0.8],'linestyle','none')
hold on;
stairs(xData(tau1:end),(mean(U_hat(burn+1:M,tau1:end))>0.5),'r-');%ylim([-5 5]);
hold off;
ylim([-0.01,1.1]);
ax = gca;
ax.XTick = xData;
xticks((xData(tau1):800:xData(end)));
datetick('x','dd/mm/yy','keepticks')
legend('$P(U_t=1|y_{1:T})$','$\hat{U}_t$','interpreter','latex')
set(gca,'fontsize',12);
print(gcf,'VolHat2.eps','-depsc');
%%
figure(123)
yyaxis left;
tt=floor(T/2);
sl=(mean(U_hat(burn+1:M,tau1:tt))>0.5);
area(xData(tau1:tt),mean(U_hat(burn+1:M,tau1:tt)),'facecolor',[0.85,0.8,0.8],'linestyle','none')
hold on;
stairs(xData(tau1:tt),sl,'r-');
ylim([-0.1,3.1]);
%
yyaxis right;
rgT=(tau1:tt);
rgTsl=rgT(sl);
plot(xData(rgT),abs(X(rgT)),'color',[0 0 0]);
hold on;
plot(xData(rgTsl),abs(X(rgTsl)),'color',[1 0 0],'LineWidth',0.1,'Marker','.','MarkerSize',10,'LineStyle','none');
ylim([-2.01 3]);
hold off;
ax = gca;
ax.XTick = xData;
xticks((xData(tau1):400:xData(tt)));
datetick('x','dd/mm/yy','keepticks')
legend('$P(U_t=1|y_{1:T})$','$\hat{U}_t$','$|y_t|$','interpreter','latex')
set(gca,'fontsize',10);
print(gcf,'VolHat2a.eps','-depsc');
%%
figure(124)
yyaxis left;
tt=floor(T/2);
sl=(mean(U_hat(burn+1:M,tt+1:end))>0.5);
area(xData(tt+1:end),mean(U_hat(burn+1:M,tt+1:end)),'facecolor',[0.85,0.8,0.8],'linestyle','none')
hold on;
stairs(xData(tt+1:end),sl,'r-');%ylim([-5 5]);
ylim([-0.1,3.1]);
%
yyaxis right;
rgT=(tt+1:T);
rgTsl=rgT(sl);
plot(xData(rgT),abs(X(rgT)),'color',[0 0 0]);
hold on;
plot(xData(rgTsl),abs(X(rgTsl)),'color',[1 0 0],'LineWidth',0.1,'Marker','.','MarkerSize',10,'LineStyle','none');
ylim([-2.01 3]);
hold off;
ax = gca;
ax.XTick = xData;
xticks((xData(tt+1):400:xData(end)));
datetick('x','dd/mm/yy','keepticks')
legend('$P(U_t=1|y_{1:T})$','$\hat{U}_t$','$|y_t|$','interpreter','latex')
set(gca,'fontsize',10);
print(gcf,'VolHat2b.eps','-depsc');
%%
scatter(alf0_hat(burn+1:end),sig20_hat(burn+1:end),'Marker','.','MarkerFaceColor','red');
xlabel('$\alpha_j$','interpreter','latex');
ylabel('$\sigma^2_j$','interpreter','latex');
hold on;
scatter(alf1_hat(burn+1:end),sig21_hat(burn+1:end),'Marker','.','MarkerFaceColor','blue');
scatter(mean(alf0_hat(burn+1:end)),mean(sig20_hat(burn+1:end)),'Marker','o','MarkerEdgeColor','black','MarkerFaceColor','black');
scatter(mean(alf1_hat(burn+1:end)),mean(sig21_hat(burn+1:end)),'Marker','o','MarkerEdgeColor','black','MarkerFaceColor','black');
hold off;
legend('Regime 0', 'Regime 1');
set(gca,'fontsize',12);
print(gcf,'Scatter.eps','-depsc');
%%
disp('Posterior mean (i.e., Bayesian estimates)');
disp([mean(alf0_hat(burn+1:end,1)),mean(alf1_hat(burn+1:end,1)),mean(sig20_hat(burn+1:end,1)),mean(sig21_hat(burn+1:end,1)),mean(theta_hat(burn+1:end,1))])
disp([mean(beta0_hat(burn+1:end,1)), mean(beta1_hat(burn+1:end,1))])