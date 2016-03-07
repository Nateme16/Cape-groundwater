%% This is the start of the groundwater transit time problem
clear all

%parameters
n=30                     % years
s=3                    % sources
k=.02                     % annual attenuation rate 
r=.1                    % discount rate


R=    .3            %net annual removel rate
IN0=100
%tau=randi([0 100],[s,1]) % distribution of transit times
tau= [3, 2, 1]'
taum=max(tau)
c=zeros(s,n+taum)         % matrix of sources

c(1,:)=10                % yearly N input rate
c(2,:)=10
c(3,:)=10

cost_s=1000
cost_in=1000
%balance=.5
control_s=ones(s+1,n)
control_s(:,:)=.1

[npv damage IN net] = estuaryi2(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0);

mnpv=@(control_s) estuaryi2(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0);

%options=optimset('Display','final','PlotFcns',@optimplotfval,'MaxFunEvals',10e10,'Algorithm','sqp');
lb=zeros(size(control_s))
ub=zeros(size(control_s))
ub(:,:)=1
%[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);
%V

options=optimset('Display','final','Algorithm','sqp','MaxFunEvals', 10e5,'TolX',10e-5);
lb=zeros(size(control_s))
ub=ones(size(control_s))

[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);

[npvop damageop INop netop]=estuaryi2(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0)
;

plot(x(:,1+taum:end)')

%plot(damage(1+taum:end))
