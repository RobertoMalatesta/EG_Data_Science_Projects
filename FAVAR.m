clear
clc

%This program takes 30 data series from the Federal Reserve Bank of St. Louis
%and implements a factor-augmented vector autogressions (FAVAR) in order
%to forecast GDP using a change-of-basis dimension reduction technique,
%a la Bernanke & Boivin (2005).

load faVAR.mat

Y_=data';


periods=[1, 2:144];

figure
plot(periods,Y_)
xlim([periods(1) periods(end)])
title('U.S. Raw Data')

%We make our data stationary by transforming them into log growth rates
%of the original series
dFedFunds=Y_(1,2:144);


lGDP=log(Y_(2,:));
dGDP=diff(lGDP);

DEFLATOR=Y_(3,:);
lDEFLATOR=log(DEFLATOR);
inflation=diff(lDEFLATOR);

dAAA=Y_(4,2:144);

lMonBase=log(Y_(5,:));
dMonBase=diff(lMonBase);

dCAbal=Y_(6,2:144);

lFxCapCon=log(Y_(7,:));
dFxCapCon=diff(lFxCapCon);

dTRfyr=Y_(8,2:144);

dTRtyr=Y_(9,2:144);

lDispInc=log(Y_(10,:));
dDispInc=diff(lDispInc);

dExRaNa=Y_(11,2:144);

dCivEmp=Y_(12,2:144);

lExRes=log(Y_(13,:));
dExRes=diff(lExRes);

lGovCon=log(Y_(14,:));
dGovCon=diff(lGovCon);

lGDPpot=log(Y_(15,:));
dGDPpot=diff(lGDPpot);

lGNP=log(Y_(16,:));
dGNP=diff(lGNP);

lGDPrInv=log(Y_(17,:));
dGDPrInv=diff(lGDPrInv);

dTRoyr=Y_(18,2:144);

dMRTGthryr=Y_(19,2:144);

dBkPr=Y_(20,2:144);

lTotEmp=log(Y_(21,:));
dTotEmp=diff(lTotEmp);

lExRaBr=log(Y_(22,:));
dExRaBr=diff(lExRaBr);

lPerSav=log(Y_(23,:));
dPerSav=diff(lPerSav);

lPerCon=log(Y_(24,:));
dPerCon=diff(lPerCon);

lPPI=log(Y_(25,:));
dPPI=diff(lPPI);

dUnEmpDur=Y_(26,2:144);

dUnEmpR=Y_(27,2:144);

lDurGCon=log(Y_(28,:));
dDurGCon=diff(lDurGCon);

lNdurGCon=log(Y_(29,:));
dNdurGCon=diff(lNdurGCon);

lMone=log(Y_(30,:));
dMone=diff(lMone);

%Finally, order the matrix
Y = [dFedFunds; dGDP; inflation; dAAA; dMonBase; dCAbal; dFxCapCon; dTRfyr;
	dTRtyr; dDispInc; dExRaNa; dCivEmp; dExRes; dGovCon; dGDPpot; dGNP;
	dGDPrInv; dTRoyr; dMRTGthryr; dBkPr; dTotEmp; dExRaBr; dPerSav;
	dPerCon; dPPI; dUnEmpDur; dUnEmpR; dDurGCon; dNdurGCon; dMone];


[N,T] = size(Y);



%Transform data by computing Z-scores

Means = zeros(30);
for i=1:30;
   Means(i)=mean(Y(i,:));
end;

Stdevs = zeros(30);
for i=1:30;
   Stdevs(i)=std(Y(i,:));
end;

Ystandard = zeros(N,T);

for i=1:30;
   for k=1:T-1;
       Ystandard(i,k) = (Y(i,k)-Means(i))/Stdevs(i);
   end;
end;

Dates = [1, 2:143];

figure
set(gcf,'Color',[1 1 1])
plot(Dates,Ystandard)
title('U.S. Data','Fontsize',14)
xlim([Dates(1) Dates(end)])
set(gca,'XTick',Dates(1):8:Dates(end))

%In order to drop the last 2 years for our model tests

Y_1 = Ystandard(:,1:T-8);



%Make a Scree plot to find the proportion of variance explained by
%eigenvalues

[W eigvals] = eig((1/(N*(T-8)))*Y_1*Y_1');
lambda=diag(eigvals);


for i=1:length(lambda);
    S(i)=sum(lambda(1:i));
end

figure
subplot(2,1,1);
plot(lambda', 'g', 'linewidth', 2);
title('Scree Plot for U.S. Data')
subplot(2,1,2);
plot(S','b', 'linewidth', 2);


%Find principal components

F = W'*Y_1;


%Estimate a VAR(1)

k=5;

PC = F(1:k,:);

LPC=PC(:,1:T-9);

A = (PC(:,2:T-8)*LPC')/(LPC*LPC');

Sigma_U=(1/(T-9))*PC(:,(2:T-8))*(eye(T-9)-(LPC'*(inv(LPC*LPC')*LPC)))*(PC(:,(2:T-8))');




%Plot the (in-sample) predictions of GDP growth, the bank rate, inflation

%One factor

 Yhat = W(:,1)*PC(1,:);
 for i=1:3;
 figure
 plot(1:T-8,Y_1(i,:),'b',1:T-8,Yhat(i,:),'g', 'linewidth', 2)
 title('In-sample 1 Factor FAVAR')
 end
 legend('Actuals','Forecast')

%Three factors
 Yhat = W(:,1:3)*PC(1:3,:);
 for i=1:3;
 figure
 plot(1:T-8,Y_1(i,:),'b',1:T-8,Yhat(i,:),'g', 'linewidth', 2)
 title('In-sample 3 Factor FAVAR')
 end
 legend('Actuals','Forecast')

%Five factors
 Yhat = W(:,1:5)*PC(1:5,:);
for i=1:3;
figure
 plot(1:T-8,Y_1(i,:),'b',1:T-8,Yhat(i,:),'g', 'linewidth', 2)
 title('In-sample 5 Factor FAVAR')
end
legend('Actuals','Forecast')


%Make (out-of-sample) forecasts

PCf = zeros(5,8);
Yf = zeros(3,8);
 
PCf(:,1) = PC(:,T-8);

for i=1:8;
PCf(:,i+1) = A*PCf(:,i);
Yf(:,i+1) = W(3,1:5)*PCf(:,i+1);
end;


Means= zeros(3);
for i=1:3;
    Means(i)=mean(Y(i,1:T-8));
end;

Stdevs = zeros(3);
for i=1:3;
    Stdevsold(i)=std(Y(i,1:T-8));
end;

for i=1:3;
    for k=1:8;
        Yf(i,k) = Means(i)+(Stdevs(i)*Yf(i,k));
    end;
end;

figure
plot(1:8,Yf(1,1:8),'g',1:8,Y(1,T-7:T),'b', 'linewidth', 2);
title('Forecasted Federal Funds Rate vs. Actual Federal Funds Rate');
legend('Forecast','Actuals')


figure
plot(1:8,Yf(2,1:8),'g',1:8,Y(1,T-7:T),'b', 'linewidth', 2);
title('Forecasted GDP Growth vs. Actual GDP Growth');
legend('Forecast','Actuals')


figure
plot(1:8,Yf(3,1:8),'g',1:8,Y(3,T-7:T),'b');
title('Forecasted Inflation vs. Actual Inflation', 'linewidth', 2);
legend('Forecast','Actuals')





%Forecast Comparisons



% Definition of univariate time series for GDP
GDP=Y(2,:);

% Define X matrix with a constant in the AR(1) process
X=[ones(1,T-1); GDP(:,1:T-1)];

% Estimate by OLS
B=(GDP(:,2:T)*X')/(X*X');
v=B(:,1);
Rho=B(:,2);

% Forecast
GDPforecast(1,1) = GDP(1,T);
for i=1:8;
    GDPforecast(1,i+1)=v+Rho*GDPforecast(1,i);
end;

% Undo transformation
for i=1:9;
    GDPforecast(1,i)=Means(3)+(Stdevs(3)*GDPforecast(1,i));
end;

clear X;
clear B;
clear v;
clear Rho;


CPI=Y(3,:);

X=[ones(1,T-1); CPI(:,1:T-1)];


B=(CPI(:,2:T)*X')/(X*X');
v=B(:,1);
Rho=B(:,2);


CPIforecast(1,1) = CPI(1,T);
for i=1:8;
    CPIforecast(1,i+1)=v+Rho*CPIforecast(1,i);
end;

for i=1:9;
    CPIforecast(1,i)=Means(2)+(Stdevs(2)*CPIforecast(1,i));
end;

clear X;
clear B;
clear v;
clear Rho;


FFR=Y(1,:);

X=[ones(1,T-1); FFR(:,1:T-1)];

% Estimation
B=(FFR(:,2:T)*X')/(X*X');
v=B(:,1);
Rho=B(:,2);


FFRforecast(1,1) = FFR(1,T);
for i=1:8;
    FFRforecast(1,i+1)=v+Rho*FFRforecast(1,i);
end;


for i=1:9;
    FFRforecast(1,i)=Means(1)+(Stdevs(1)*FFRforecast(1,i));
end;

clear X;
clear B;
clear v;
clear Rho;



%Predicted values, trivariate VAR(p)

VAR(1,:)=Y(1,:);
VAR(2,:)=Y(2,:);
VAR(3,:)=Y(3,:);
p=2;
y=VAR(:,p+1:end);


for j=1:p
    X=[ones(1,T-p); VAR(:,p+1-j:end-j);];
end

B=y*X'/(X*X');
v=B(:,1);
Rho=B(:,2:4);

% Forecast
VARforecast(:,1)=VAR(:,T);
for i=1:8;
    VARforecast(:,i+1)=v+Rho*VARforecast(:,i);
end;


for i=1:9;
    VARforecast(1,i)=Means(1)+(Stdevs(1)*VARforecast(1,i));
    VARforecast(2,i)=Means(2)+(Stdevs(2)*VARforecast(2,i));
    VARforecast(3,i)=Means(3)+(Stdevs(3)*VARforecast(3,i));
end;


figure
subplot(3,1,1);
plot(1:9,Yf(1,1:9),'g',1:9,Y(1,T-8:T),'b',1:9,FFRforecast,'r', 'linewidth', 3);
title('FedFunds Rate vs. Forecasts');
legend('Factor Model','Actuals','AR(1)','Location','Southwest');
subplot(3,1,2);
plot(1:9,Yf(2,1:9),'g',1:9,Y(2,T-8:T),'b',1:9,GDPforecast,'r', 'linewidth', 3);
title('Inflation vs. Forecasts');
legend('Factor Model','Actuals','AR(1)','Location','Southwest');
subplot(3,1,3);
plot(1:9,Yf(3,1:9),'g',1:9,Y(3,T-8:T),'b',1:9,CPIforecast,'r', 'linewidth', 3);
title('GDP Growth Rate vs. Forecasts');
legend('Factor Model','Actuals','AR(1)','Location','Southwest');



figure
subplot(3,1,1);
plot(1:9,Yf(1,1:9),'g',1:9,Y(1,T-8:T),'b',1:9,VARforecast(1,:),'r', 'linewidth', 3);
title('FedFunds Rate vs. Forecasts');
legend('Factor Model','Actuals','VAR(2)','Location','Southwest');
subplot(3,1,2);
plot(1:9,Yf(2,1:9),'g',1:9,Y(2,T-8:T),'b',1:9,VARforecast(2,:),'r', 'linewidth', 3);
title('Inflation vs. Forecasts');
legend('Factor Model','Actuals','VAR(2)','Location','Southwest');
subplot(3,1,3);
plot(1:9,Yf(3,1:9),'g',1:9,Y(3,T-8:T),'b',1:9,VARforecast(3,:),'r', 'linewidth', 3);
title('GDP Growth Rate vs. Forecasts');
legend('Factor Model','Actuals','VAR(2)','Location','Southwest');