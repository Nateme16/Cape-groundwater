
%% This is the start of the groundwater transit time problem

%Author: Nate Merrill

%this file has the input paramters and calls the function files estuaryi2
%which is the model of the estuary through time, which uses a file dettman
%(name of a scientist (Ed Dettman) here whose "model" I adapted  for estuary N stocks)

clear all

%parameters- these are all placeholders until we parameterize it for an
%estuary system

n=40                   % years
s=4                    % sources
k=0                  % annual attenuation rate 
r=.05                  % discount rate
B=  20983 +21561       % background/ non controllable sources
R=   0                % net N caryover from year before
Target= 52000          % target N load
gamma=100                % penalty for difference
shellmax= 1400              % max N removal from aquiculture
IN0=1                 %initial N stock?

%tau=randi([0 100],[s,1]) % distribution of transit times
tau= [1,5,10,15,20]'        %transit time for each source
taum=max(tau)          %this is used to deal with load before choice years
c=zeros(s,n+taum)      % matrix of sources

c(1,:)=5427% yearly N input rate
c(2,:)=7562
c(3,:)=6885
c(4,:)=5188
c(5,:)=2176



cost_s=31             % cost parameter for source control
cost_in=10           % cost parameter for in-estuary control

control_s=zeros(s+1,n);  % each cell represents the total load abated at source
     
% This part gets the optimal control back and other outputs
[npv damage IN net cost kgs kgi] = estuaryi4(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma);

mnpv=@(control_s) estuaryi4(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma);

options=optimset('Display','final','Algorithm','sqp','MaxFunEvals', 10e5,'TolX',10e-10,'PlotFcns',@optimplotfval,'MaxIter',1000);
%options=optimset('Display','final','Algorithm','sqp','PlotFcns',@optimplotfval);

lb=zeros(size(control_s));
ub=ones(size(control_s));

for i=1:s;
   ub(i,:)= c(i,1);
end
ub(s+1,:)=shellmax

[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,[],options);

[npvop damageop INop netop cost kgsop kgiop]=estuaryi4(x,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma);


%plots control and in estuary N over time
subplot(2,1,1)
plot(x(:,1:end)')
subplot(2,1,2)
plot(INop(:,1+taum:end)-Target)


%plot(damage(1+taum:end))
