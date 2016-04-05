%% This is the start of the groundwater transit time problem
clear all

%parameters
n=100                     % years
s=10                     % sources
k=.02                     % annual attenuation rate 
r=.05                    % discount rate

%tau=randi([0 100],[s,1]) % distribution of transit times
tau= [50,30,20,10,5,4,2,1,1,1]'
taum=max(tau)
c=zeros(s,n+taum)         % matrix of sources

c(:,:)=10
%c(1,:)=16                % yearly N input rate
%c(2,:)=16
%c(3,:)=16

cost_s=500
%cost_in=.2
%balance=.5
control_s=ones(s,1)
control_s(:,:)=.8

[npv damage TC net ] = estuarystatic(control_s,c,taum,n,k,tau,s,r,cost_s);

mnpv=@(control_s) estuarystatic(control_s,c,taum,n,k,tau,s,r,cost_s);

options=optimset('Display','final','PlotFcns',@optimplotfval,'Algorithm','sqp');
lb=zeros(size(control_s))
ub=ones(size(control_s))

[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);
%options=optimset('Display','final','MaxIter',10e100,'MaxFunEvals',10e100,'TolFun',10e-100,'TolX',10e-100);

[npvop damageop TCop netop ]=estuarystatic(x,c,taum,n,k,tau,s,r,cost_s);

x

plot(damageop(1+taum:end))
