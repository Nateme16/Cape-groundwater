
%% This is the start of the groundwater transit time problem

%Author: Nate Merrill

%this file has the input paramters and calls the function files estuaryi2
%which is the model of the estuary through time, which uses a file dettman
%(name of a scientist (Ed Dettman) here whose "model" I adapted  for estuary N stocks)

clear all

%parameters- these are all placeholders until we parameterize it for an
%estuary system

n=20                   % years
s=3                    % sources
k=.1                   % annual attenuation rate 
r=.05                  % discount rate
B=100                  % background/ non controllable sources
R=   .5                % net N caryover from year before
IN0=100
%tau=randi([0 100],[s,1]) % distribution of transit times
tau= [5, 3, 1]'        %transit time for each source
taum=max(tau)          %this is used to deal with load before choice years
c=zeros(s,n+taum)      % matrix of sources

c(1,:)=10                % yearly N input rate
c(2,:)=10
c(3,:)=10

cost_s=1000             % cost parameter for source control
cost_in=2000            % cost parameter for in-estuary control

control_s=ones(s+1,n);  % each cell represents the % of total load left after treatment of each source
control_s(:,:)=.1;      % each cell is the in % in estuary N left in estuary after treatment

% This part gets the optimal control back and other outputs
[npv damage IN net] = estuaryi2(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B);

mnpv=@(control_s) estuaryi2(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B);

options=optimset('Display','final','Algorithm','sqp','MaxFunEvals', 10e5,'TolX',10e-5,'PlotFcns',@optimplotfval);
lb=zeros(size(control_s));
ub=ones(size(control_s));

[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);

[npvop damageop INop netop]=estuaryi2(x,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B);


%plots control and in estuary N over time
subplot(2,1,1)
plot(x(:,1+taum:end)')
subplot(2,1,2)
plot(INop(:,1+taum:end))


%plot(damage(1+taum:end))
