clear all; %Ex2.1 2017-2017
clc;
dataAll=xlsread('Cyber data 2017.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2017=dataAll(:,5);
%
T=size(data_2017,1); % number of daily observations
n=size(data_2017,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2017');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2017);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2017(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);


% total cyber attack for 2018 
clear all; %cybercrime 2017-2018
clc;
dataAll=xlsread('Cyber data 2018.xlsx');
% column 1: Cyber Crime
% column 2: Cyber Espionage
% column 3: Cyber Warfare
% column 4: Hacktivism
% column 5: Total 
%
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2018=dataAll(:,5);
%
T=size(data_2018,1); % number of daily observations
n=size(data_2018,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2018');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2018(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2018);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2018(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2018(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cyber crime 2017-2018
clear all; %cybercrime 2017-2018
clc;
dataAll=xlsread('CyberData.xlsx');
%
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data=dataAll(:,1);
%
T=size(data,1); % number of daily observations
n=size(data,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%% cyber crime 2017 

clear all; 
clc;
dataAll=xlsread('Cyber data 2017.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2017=dataAll(:,2);
%
T=size(data_2017,1); % number of daily observations
n=size(data_2017,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2017');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2017);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2017(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%%%% Cyber crime 2018 

clear all; 
clc;
dataAll=xlsread('Cyber data 2018.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2018=dataAll(:,1);
%
T=size(data_2018,1); % number of daily observations
n=size(data_2018,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2018');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2018(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2018);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2018(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2018(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hacktivism 2017-2018
clear all; %Hacktivism 2017-2018
clc;
dataAll=xlsread('CyberData.xlsx');
%
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data=dataAll(:,4);
%
T=size(data,1); % number of daily observations
n=size(data,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%% Hacktivism 2017 

clear all; 
clc;
dataAll=xlsread('Cyber data 2017.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2017=dataAll(:,2);
%
T=size(data_2017,1); % number of daily observations
n=size(data_2017,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2017');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2017);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2017(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%%%% Hacktivism 2018 

clear all; 
clc;
dataAll=xlsread('Cyber data 2018.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2018=dataAll(:,1);
%
T=size(data_2018,1); % number of daily observations
n=size(data_2018,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2018');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2018(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2018);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2018(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2018(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);



%Cyber espionage  2017-2018
clear all; %cyberespionage 2017-2018
clc;
dataAll=xlsread('CyberData.xlsx');
%
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data=dataAll(:,2);
%
T=size(data,1); % number of daily observations
n=size(data,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%% cyber espionage  2017 

clear all; 
clc;
dataAll=xlsread('Cyber data 2017.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2017=dataAll(:,2);
%
T=size(data_2017,1); % number of daily observations
n=size(data_2017,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2017');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2017);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2017(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%%%% Cyber espionage 2018 

clear all; 
clc;
dataAll=xlsread('Cyber data 2018.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2018=dataAll(:,1);
%
T=size(data_2018,1); % number of daily observations
n=size(data_2018,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2018');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2018);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2018(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2018(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);








%Cyber warfare 2017-2018
clear all; %cyberwarfare 2017-2018
clc;
dataAll=xlsread('CyberData.xlsx');
%
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data=dataAll(:,1);
%
T=size(data,1); % number of daily observations
n=size(data,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%% cyber cyberwarfare 2017 

clear all; 
clc;
dataAll=xlsread('Cyber data 2017.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2017=dataAll(:,2);
%
T=size(data_2017,1); % number of daily observations
n=size(data_2017,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2017');
endDate=datenum('12-31-2017');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2017);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2017(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2017(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

%%%% Cyberwarfare 2018 

clear all; 
clc;
dataAll=xlsread('Cyber data 2018.xlsx');
%%%%%%%%%%%%%%%%%%%%%% Select the 5th time series
data_2018=dataAll(:,1);
%
T=size(data_2018,1); % number of daily observations
n=size(data_2018,2); % number of series
%%%%%%%%%%%%%%%%%%%%%% Plot Data
startDate = datenum('01-01-2018');
endDate=datenum('12-31-2018');
xData = linspace(startDate,endDate,T);
%
figure(1)
stairs(xData,data_2017(:,1)); %plot the time series
xlim([xData(1),xData(end)]);
ax = gca;
ax.XTick = xData;
xticks((xData(1):150:xData(end)));
datetick('x','dd/mm/yy','keepticks')
set(gca,'fontsize',12);
print(gcf,'Cyber.eps','-depsc');
%%
%%%%%%%%%%%%%%%%%%%%%% Bayesian Inference
% Ga(alfa,beta) prior hyperparameters
alfa=2; 
beta=2;
dlt=0.01;
grd=(0:dlt:5);
ngrd=size(grd,2);
pdfPrior=pdf('gamma',grd,alfa,1/beta); %prior pdf
%
% Ga(alfabar,betabar) posterior parameters
alfabar=alfa+sum(data_2018);                       
betabar=beta+T;         
pdfPosterior=pdf('gamma',grd,alfabar,1/betabar); 
%
% Ga(alfabar,betabar) posterior parameters SUB-sample
alfabar1=alfa+sum(data_2018(1:30,1));                       
betabar1=beta+30;   
pdfPosterior1=pdf('gamma',grd,alfabar1,1/betabar1);
%
%%%%%%%%%%%%% Plot prior and posterior distributions
figure(2)
plot(grd,pdfPrior*dlt,'r-');
hold on;
plot(grd,pdfPosterior1*dlt,'k--');
plot(grd,pdfPosterior*dlt,'k-');
hold off;
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Prior','Posterior 30 Days','Posterior 365 Days');
set(gca,'fontsize',12);
print(gcf,'Posterior.eps','-depsc');
%
% Sequential estimation incraesing the sample size
alfabarCum=zeros(T,1);
betabarCum=zeros(T,1);
alfabarCum(1,1)=alfa;                       
betabarCum(1,1)=beta; 
for t=2:T
    alfabarCum(t,1)=alfabarCum(t-1,1)+data_2018(t,1);                       
    betabarCum(t,1)=betabarCum(t-1,1)+1;
end
%
figure(3)
plot((1:T)',alfabarCum./betabarCum,'k-');
hold on;
plot((1:T)',alfa/beta*ones(T,1),'r:');
plot((1:T)',alfabar/betabar*ones(T,1),'r--');
hold off;
xlim([1,T]);
ylim([0,5]);
legend('Sequential Estimate','Prior Mean','Posterior Mean Whole Sample');
set(gca,'fontsize',12);
print(gcf,'PosteriorSeq.eps','-depsc');
%%
lambdahat=alfabar/betabar;
q0025=icdf('gamma',0.025,alfabar,1/betabar);
q0975=icdf('gamma',0.975,alfabar,1/betabar);
sl=boolean((cumsum(pdfPosterior)*dlt>0.025).*(cumsum(pdfPosterior)*dlt<0.975));
Zoom=boolean((cumsum(pdfPosterior)*dlt>0.0001).*(cumsum(pdfPosterior)*dlt<0.9999));
%
figure(4)
plot(grd(Zoom),pdfPosterior(Zoom)*dlt,'k-');
hold on;
area(grd(sl),pdfPosterior(sl)*dlt,'facecolor',[0.9,0.9,0.9],'linestyle','none')
plot(alfabar/betabar,0,'ro');
str='$\hat{\lambda}_{B}$';
text(alfabar/betabar,-0.005,str,'interpreter','latex');
hold off;
xticks([q0025,q0975]);
set(gca,'xticklabels',{['q_{0.975}=',num2str(q0025,4)],['q_{0.975}=',num2str(q0975,4)]});
ylim([0,max(pdfPosterior*dlt*1.5)]);
legend('Posterior distribution','High Probability Density (HPD) region','Posterior mean');
set(gca,'fontsize',12);
print(gcf,'PosteriorHPD.eps','-depsc');
%%
%%%%%%%%% Output
disp(['Prior mean=',num2str(alfa/beta,4)]);
disp(['Intensity estimate (posterior mean)=',num2str(lambdahat,4)]);
disp(['Credible Interval=(',num2str(q0025,4),',',num2str(q0975,4),')']);

