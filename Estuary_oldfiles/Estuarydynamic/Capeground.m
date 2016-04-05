%% This is the start of the groundwater transit time problem
clear all

%parameters
n=20                     % years
s=3                     % sources
k=.02                     % annual attenuation rate 
r=.1                    % discount rate

%tau=randi([0 100],[s,1]) % distribution of transit times
tau= [3, 2, 1]'
taum=max(tau)
c=zeros(s,n+taum)         % matrix of sources

c(1,:)=10                % yearly N input rate
c(2,:)=10
c(3,:)=10

cost_s=200
%cost_in=.2
%balance=.5
control_s=ones(s,n)
control_s(:,:)=.8

[npv damage TC net] = estuary(control_s,c,taum,n,k,tau,s,r,cost_s);

mnpv=@(control_s) estuary(control_s,c,taum,n,k,tau,s,r,cost_s);

%options=optimset('Display','final','PlotFcns',@optimplotfval,'MaxFunEvals',10e10,'Algorithm','sqp');
lb=zeros(size(control_s))
ub=zeros(size(control_s))
ub(:,:)=1
%[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);

options=optimset('Display','final','PlotFcns',@optimplotfval,'Algorithm','sqp','MaxFunEvals', 10e10);
lb=zeros(size(control_s))
ub=ones(size(control_s))

[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);

[npvop damageop TCop netop]=estuary(x,c,taum,n,k,tau,s,r,cost_s);

plot(x(:,1+taum:end)')

%plot(damage(1+taum:end))
