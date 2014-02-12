%This program estimates the parameters of a Bayesian structural vector
%autoregression by using a Random Walk Metropolis-Hastings algortihm 
%with uniform truncated priors. 
%Y_t = \phi(Y_t-1) + e_t


clear all
clc


load Bayes.mat
Y = Bayes;


r = Y(:,2);
cpi = Y(:,3);
gdp = Y(:,1);

dates = 1955:0.25:2011;


% We first difference the logs of GDP and Inflation and then demean the data

lcpi = (log(cpi(2:end)) - log(cpi(1:end-1)))*100;
lgdp = (log(gdp(2:end)) - log(gdp(1:end-1)))*100;
r = r(2:end);
t = size(r);


mcpi = (lcpi - repmat(mean(lcpi),t,1));
mgdp = (lgdp - repmat(mean(lgdp),t,1));
mr = (r - repmat(mean(r),t,1));


figure
plot(dates(1:end-1),mcpi,'r')
hold on
plot(dates(1:end-1),mgdp)
plot(dates(1:end-1),mr,'g')
legend('CPI','GDP','mr')


% We first need to estimate using OLS to get the estimates for \phi

z = [mr mgdp mcpi];
Y = z';
[K,T] = size(Y);

P = 1; %"optimal lag length"
X = Y(:,1:T-P);
phi = Y(:,P+1:T)*X'/(X*X');
%A0 = [2 1 3; 2 8 5; 3 9 1];
%A1 = [9 7 2; 6 9 1; 1 7 2];
%phi = A0\A1;

Su=(1/(T-K*P-1))*Y(:,P+1:T)*(eye(T-P)-X'*inv(X*X')*X)*Y(:,P+1:T)';

%These are the initial values for the vector of parameters \theta

theta =[phi(1,1:3) phi(2,1:3) phi(3,1:3) Su(:,1)' Su(2:3,2)' Su(3,3)]';

LB = -2*ones(size(theta));
UB = 2*ones(size(theta));
x = theta;

MSE = zeros(K,K,T);
MSE(:,:,1)=Su; 

for j=2:T
    MSE(:,:,j)=MSE(:,:,j-1)+(phi^(j-1))*Su*(phi^(j-1))';
end
sigma_yy = MSE(:,:,224);

% Number of draws and then we initialize the proposal variance

J = 2000;
epseye=0.1;

bdraw=x';
lpostdraw = -9e+200;
vscale=diag(abs(theta))*epseye+1e-4*eye(length(x)); 

%The above is the variance covariance matrix of the proposal distribution

bb_=zeros(length(x),J);

%These matrices that keep track of switches and draws outside LB and UB
OutsideProp=zeros(J,1); 
SwitchesProp=zeros(J,1);

%Number of draws that do not satisfy the stability of phi and boundary conditions
q=0;

%Number of switches (acceptances)
pswitch=0;

%Iteration counter
iter=0;


% Implementation of the Metropolis-Hastings algorithm

tic
for iter=1:J
    
    if mod(iter,0.05*J)==0
        disp(iter)
        toc
    end
    
    
    iter=iter+1;
    %The draw is from proposal density Theta*_{t+1} ~ N(Theta_{t},vscale)
    
    bcan = mvnrnd(bdraw, vscale,1)'; %where 1 is an indicator for 1 draw - draw vectors separate
    newdr = [bcan(10) 0 0 ; bcan(11) bcan(13) 0; bcan(12) bcan(14) bcan(15)];
    
    newdrp = [bcan(1:3)'; bcan(4:6)'; bcan(7:9)'];
    if max(abs(eig(newdrp)))<1
          
    if min(bcan > LB)==1 %need to be inside the interval
        if min(bcan < UB)==1 %need to be inside the interval
            %lpostcan = LLDSGE(bcan,Z);
            lpostcan = LOGPRIOR_U_f(bcan)+LLVAR(bcan,X,z);%switch on for use of priors
            %lpostcan = LOGPRIOR_U_f(bcan);%switch on for prior predictive analysis
            laccprob = lpostcan-lpostdraw; %compare the loglikelihood with the current draw
        else
            laccprob=-9e+200;
            q=q+1;
        end
    else
        laccprob=-9e+200;
        q=q+1;
    end
    else
        laccprob=-9e+200;
        q = q+1;
    end
    
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        bdraw=bcan;
        pswitch=pswitch+1; %if accepted, then 1 draw
    end
    
    
    bb_(:,iter)=bcan; %store the draw (remember the bb_ store vector)
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    
    
end
toc

disp(['iter: ',num2str(iter)]); %number of acceptances
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]); %number of rejections

figure
bb_=bb_(:,2:end);
for j=1:15;
    subplot(3,5,j); 
    hist(bb_(j,:),50);
end

convcheck(bb_);

figure
plotpost(bb_,0) %=0 then we plot the pdf instead of the cumulative, see plotpo
%% Probability intervals

%Bayesian Estimation of the Unrestricted SVAR

clear all
clc

load Bayes.mat
Y = Bayes;


r = Y(:,2);
cpi = Y(:,3);
gdp = Y(:,1);

dates = 1955:0.25:2011;

lcpi = (log(cpi(2:end)) - log(cpi(1:end-1)))*100;
lgdp = (log(gdp(2:end)) - log(gdp(1:end-1)))*100;
r = r(2:end);
t = size(r);


mcpi = (lcpi - repmat(mean(lcpi),t,1));
mgdp = (lgdp - repmat(mean(lgdp),t,1));
mr = (r - repmat(mean(r),t,1));


figure 
plot(dates(1:end-1),mcpi,'r')
hold on
plot(dates(1:end-1),mgdp)
plot(dates(1:end-1),mr,'g')
legend('CPI','GDP','mr')


z = [mr mgdp mcpi];
Y = z';
[K,T] = size(Y);

P = 1; 
X = Y(:,1:T-P);

%such that phi = Y(:,P+1:T)*X'/(X*X');

A0 = randn(3);
A1 = randn(3);
A0 = A0*0.7;
A1 = A1*0.05;
phi = A0\A1;

%Su=(1/(T-K*P-1))*Y(:,P+1:T)*(eye(T-P)-X'*inv(X*X')*X)*Y(:,P+1:T)';

%%
theta =[A0(1,1:3) A0(2,1:3) A0(3,1:3) A1(1,1:3) A1(2,1:3) A1(3,1:3)]'; 

LB = -2*ones(size(theta));
UB = 2*ones(size(theta));
x = theta;

%MSE = zeros(K,K,T);
%MSE(:,:,1)=Su; 

%for j=2:T
%    MSE(:,:,j)=MSE(:,:,j-1)+(phi^(j-1))*Su*(phi^(j-1))';
%end
%sigma_yy = MSE(:,:,224);

% Number of draws
J = 20000;
epseye=0.0001;

bdraw=x';
lpostdraw = -9e+200;
vscale=diag(abs(theta))*epseye+1e-4*eye(length(x)); 

bb_=zeros(length(x),J);

OutsideProp=zeros(J,1); 
SwitchesProp=zeros(J,1);

q=0;
%Number of switches (acceptances)
pswitch=0;
%Iteration counter
iter=0;

%Implementation of the M-H algorithm

tic
for iter=1:J
    
    if mod(iter,0.05*J)==0
        disp(iter)
        toc
    end
    
    iter=iter+1; 
    bcan = mvnrnd(bdraw, vscale,1)'; %where 1 is an indicator for 1 draw - draw vectors separate
    %newdre = [bcan(10) 0 0 ; bcan(11) bcan(13) 0; bcan(12) bcan(14) bcan(15)]; 
    
    a0 = [bcan(1:3)'; bcan(4:6)'; bcan(7:9)'];
    phi = [bcan(1:3)'; bcan(4:6)'; bcan(7:9)']\[bcan(10:12)'; bcan(13:15)'; bcan(16:18)'];
    Su = [bcan(1:3)'; bcan(4:6)'; bcan(7:9)']*[bcan(1:3)'; bcan(4:6)'; bcan(7:9)']';
    
    if max(abs(eig(phi)))<1
         
    if min(bcan > LB)==1 
        if min(bcan < UB)==1 
            if det(a0)~=0
            %lpostcan = LLDSGE(bcan,Z);
            lpostcan = LOGPRIOR_U(bcan)+LLVAR_df(X,z,phi,Su);
            %lpostcan = LOGPRIOR(bcan);%switch on for prior predictive analysis
            laccprob = lpostcan-lpostdraw; 
            else
                laccprob=-9e+200;
            q=q+1;
            end
        else
            laccprob=-9e+200;
            q=q+1;
        end
    else
        laccprob=-9e+200;
        q=q+1;
    end
    else
        laccprob=-9e+200;
        q=q+1;
    end
    
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        bdraw=bcan;
        pswitch=pswitch+1; %if accepted, then 1 draw
    end
    
    bb_(:,iter)=bdraw; %store the draw (remember the bb_ store vector)
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    
    
end
toc

disp(['iter: ',num2str(iter)]); %number of acceptances
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]); %number of rejections

figure
bb_=bb_(:,2:end);
for j=1:18;
    subplot(3,6,j); %2x4 because we have 8 parameters
    hist(bb_(j,:),50);
end

convcheck(bb_);

figure
plotpost(bb_,0) %=0 then we plot the pdf instead of the cumulative, see plotpo


%Bayesian Estimation of the Restricted SVAR

clear all
clc

load Bayes.mat
Y = Bayes;


r = Y(:,2);
cpi = Y(:,3);
gdp = Y(:,1);

dates = 1955:0.25:2011;

lcpi = (log(cpi(2:end)) - log(cpi(1:end-1)))*100;
lgdp = (log(gdp(2:end)) - log(gdp(1:end-1)))*100;
r = r(2:end);
t = size(r);


mcpi = (lcpi - repmat(mean(lcpi),t,1));
mgdp = (lgdp - repmat(mean(lgdp),t,1));
mr = (r - repmat(mean(r),t,1));


figure 
plot(dates(1:end-1),mcpi,'r', 'linewidth', 2.5)
hold on
plot(dates(1:end-1),mgdp, 'linewidth', 2.5)
plot(dates(1:end-1),mr,'g', 'linewidth', 2.5)
legend('CPI','GDP','mr')



z = [mcpi mgdp mr];
Y = z';
[K,T] = size(Y);

P = 1; %"optimal lag length"
X = Y(:,1:T-P);
%phi = Y(:,P+1:T)*X'/(X*X');
%A0 = [2 0 0; 2 3 0; 3 9 1];
%A1 = [9 7 2; 6 9 1; 1 7 2];
phi = Y(:,P+1:T)*X'/(X*X');
Su=(1/(T-K*P-1))*Y(:,P+1:T)*(eye(T-P)-X'*inv(X*X')*X)*Y(:,P+1:T)';
chol = chol(Su,'lower');
A0 = inv(chol);
A1 =chol*phi;
%phi = A0\A1;

%Su=(1/(T-K*P-1))*Y(:,P+1:T)*(eye(T-P)-X'*inv(X*X')*X)*Y(:,P+1:T)';

%%
theta =[A0(1,1) A0(2,1:2) A0(3,1:3) A1(1,1:3) A1(2,1:3) A1(3,1:3)]'; 

LB = -2*ones(size(theta));
UB = 2*ones(size(theta));
x = theta;

%MSE = zeros(K,K,T);
%MSE(:,:,1)=Su; 

%for j=2:T
%    MSE(:,:,j)=MSE(:,:,j-1)+(phi^(j-1))*Su*(phi^(j-1))';
%end
%sigma_yy = MSE(:,:,224);

% Number of draws
%J = 2e3;
J = 50000;
epseye=1e-4;
%Initializes the proposal variance
bdraw=x';
lpostdraw = -9e+200;
vscale=diag(abs(theta))*epseye+1e-4*eye(length(x)); 


bb_=zeros(length(x),J);


OutsideProp=zeros(J,1); %a counter for how many draws are outside the bounds
SwitchesProp=zeros(J,1);
%Number of draws that does not satisfies our condition (stability of phi
%and boundary conditions)
q=0;
%Number of switches (acceptances)
pswitch=0;
%Iteration counter
iter=0;



tic
for iter=1:J
    
    if mod(iter,0.05*J)==0
        disp(iter)
        toc
    end
    
    iter=iter+1; 
    %bcan = bdraw + norm_rnd(vscale);
    bcan = mvnrnd(bdraw, vscale,1)'; %where 1 is an indicator for 1 draw - draw vectors separate
    %newdre = [bcan(10) 0 0 ; bcan(11) bcan(13) 0; bcan(12) bcan(14) bcan(15)]; 
    a0 = [bcan(1) 0 0; bcan(2:3)' 0; bcan(4:6)'];
    phi = [bcan(1) 0 0; bcan(2:3)' 0; bcan(4:6)']\[bcan(7:9)'; bcan(10:12)'; bcan(13:15)'];
    Su = inv([bcan(1) 0 0; bcan(2:3)' 0; bcan(4:6)'])*(inv([bcan(1) 0 0; bcan(2:3)' 0; bcan(4:6)']))';
    
    if max(abs(eig(phi)))<1
        
    if min(bcan > LB)==1 %need to be inside the interval
        if min(bcan < UB)==1 %need to be inside the interval
          if det(a0)~=0 
            %lpostcan = LLDSGE(bcan,Z);
            lpostcan = LOGPRIOR_U(bcan)+LLVAR_df(X,z,phi,Su);%switch on for use of priors
            %lpostcan = LOGPRIOR(bcan);%switch on for prior predictive analysis
            laccprob = lpostcan-lpostdraw; %compare the loglikelihood with the current draw (slide 15 in bayesinad DSGE)
        else
            laccprob=-9e+200;
            q=q+1;
          end
        else
            laccprob=-9e+200;
            q=q+1;
        end
    else
        laccprob=-9e+200;
        q=q+1;
    end
    else
        laccprob=-9e+200;
        q = q+1;
    end
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        bdraw=bcan;
        pswitch=pswitch+1; %if accepted, then 1 draw
    end
    
    bb_(:,iter)=bdraw; %store the draw (remember the bb_ store vector)
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    
    
end
toc

disp(['iter: ',num2str(iter)]); %number of acceptances
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]); %number of rejections

figure
bb_=bb_(:,2:end);
for j=1:15;
    subplot(3,5,j); %2x4 because we have 8 parameters
    hist(bb_(j,:),50);
end

convcheck(bb_);

figure
plotpost(bb_,0) 


% Plotting the Impulse Responses


g=36; %(length of IRFs);

%Calculate the predicted impulse response
for i=1:(g-1)
temp=phi^i;                     %taking the coefficient matrix to power 1
impulse(:,:,i+1)=temp*a0;
end
impulse(:,:,1)=a0;


%Now draw alternative impulse responses

k=100;

for z=1:k; %outer loop over random draws
temp1=randi(J,1); %determine which column to draw from bb_
temp2=bb_(:,temp1); %draw the random coefficient vector from bb
temp_a0= [temp2(1) 0 0; temp2(2:3)' 0; temp2(4:6)'];
temp_phi= [temp2(1) 0 0; temp2(2:3)' 0; temp2(4:6)']\[temp2(7:9)'; temp2(10:12)'; temp2(13:15)'];

for i=1:g-1
temp=temp_phi^i; 
impulse_conf(:,:,i+1,z)=temp*temp_a0; %calculate impulse reponses
end
impulse_conf(:,:,1,z)=temp_a0;
end


impulse_conf_unsorted=sort(impulse_conf,4); 
impulse_conf=sort(impulse_conf,4); 
lower=floor(0.025*z);
median=round(0.5*z);
upper=ceil(0.975*z);
imp_lower=impulse_conf(:,:,:,lower); %lower bound
imp_median=impulse_conf(:,:,:,median); %median
imp_upper=impulse_conf(:,:,:,upper); %upper bound

figure

plot(squeeze(imp_median(1,3,:)),'r');
title('Response of GDP to Interest Rate')
hold on
plot(squeeze(imp_lower(1,3,:)),'b');
plot(squeeze(imp_upper(1,3,:)),'b');

figure
plot(squeeze(imp_median(2,3,:)),'r');
title('Response of Inflation to Interest Rate')
hold on
plot(squeeze(imp_lower(2,3,:)),'b');
plot(squeeze(imp_upper(2,3,:)),'b');

figure
plot(squeeze(imp_median(3,3,:)),'r');
title('Response of Interest Rate to Interest Rate')
hold on
plot(squeeze(imp_lower(3,3,:)),'b');
plot(squeeze(imp_upper(3,3,:)),'b');


% Variance Decomposition Median


imp_sq_sum=zeros(g,9,k);
for z=1:k
imp_u1(:,:,z)=impulse_conf_unsorted(:,1,:,z);
imp_u2(:,:,z)=impulse_conf_unsorted(:,2,:,z);
imp_u3(:,:,z)=impulse_conf_unsorted(:,3,:,z);
imp_flat(:,:,z)=[imp_u1(:,:,z)' imp_u2(:,:,z)' imp_u3(:,:,z)'];
imp_sq(:,:,z)=imp_flat(:,:,z).^2;
for i=1:g
    imp_sq_sum(i+1,:,z)=imp_sq_sum(i,:,z)+imp_sq(i,:,z);
end
MSE1(:,:,z)=imp_sq_sum(:,1,z)+imp_sq_sum(:,4,z)+imp_sq_sum(:,7,z);
MSE2(:,:,z)=imp_sq_sum(:,2,z)+imp_sq_sum(:,5,z)+imp_sq_sum(:,8,z);
MSE3(:,:,z)=imp_sq_sum(:,3,z)+imp_sq_sum(:,6,z)+imp_sq_sum(:,9,z);

%variance decomposition/total variance

decomp_y1(:,1,z)=imp_sq_sum(:,1,z)./MSE1(:,:,z);
decomp_y1(:,2,z)=imp_sq_sum(:,4,z)./MSE1(:,:,z);
decomp_y1(:,3,z)=imp_sq_sum(:,7,z)./MSE1(:,:,z);
decomp_y2(:,1,z)=imp_sq_sum(:,2,z)./MSE2(:,:,z);
decomp_y2(:,2,z)=imp_sq_sum(:,5,z)./MSE2(:,:,z);
decomp_y2(:,3,z)=imp_sq_sum(:,8,z)./MSE2(:,:,z);
decomp_y3(:,1,z)=imp_sq_sum(:,3,z)./MSE3(:,:,z);
decomp_y3(:,2,z)=imp_sq_sum(:,6,z)./MSE3(:,:,z);
decomp_y3(:,3,z)=imp_sq_sum(:,9,z)./MSE3(:,:,z);
end

decomp=[decomp_y1 decomp_y2 decomp_y3];
decomp_sort=sort(decomp,3);
decomp_sort=decomp_sort(2:end,:,:);

decomp_lower=decomp_sort(end,:,lower);
decomp_median=decomp_sort(end,:,median);
decomp_upper=decomp_sort(end,:,upper);

decomp_conf=[decomp_lower;decomp_median;decomp_upper];
