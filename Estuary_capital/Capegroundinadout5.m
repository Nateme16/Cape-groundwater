
%% This is the start of the groundwater transit time problem

%Author: Nate Merrill

%this file has the input paramters and calls the function files estuaryi2
%which is the model of the estuary through time, which uses a file dettman
%(name of a scientist (Ed Dettman) here whose "model" I adapted  for estuary N stocks)

clear all

%parameters- these are all placeholders until we parameterize it for an
%estuary system

n=50                   % years
s=5                    % sources
k=.1                  % annual attenuation rate 
r=.05                % discount rate
B=  20983 +21561       % background/ non controllable sources
R=   0               % net N caryover from year before
Target= 52000          % target N load
gamma=1                % penalty for difference
shellmax= 1400         % max N removal from aquiculture
maxsept= .6             %max % removal
IN0=1                 %initial N stock?

%tau=randi([0 100],[s,1]) 
% distribution of transit times
tau= [1,5,10,15,20]'        %transit time for each source
taum=max(tau)          %this is used to deal with load before choice years
cs=zeros(s,n+taum)      % matrix of sources

cs(1,:)=5427  % yearly N input rate
cs(2,:)=7562
cs(3,:)=6885
cs(4,:)=5188
cs(5,:)=2176

cost_s=1000             % cost parameter for source control
cost_in=100          % cost parameter for in-estuary control

control_s=zeros(s+1,n);  % each cell represents the total load abated at source
control_s(:,:)=1   
% This part gets the optimal control back and other outputs
[npv damage IN net cost kgs kgi scap] = estuaryi5(control_s,cs,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma);

mnpv=@(control_s) estuaryi5(control_s,cs,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma);

options=optimset('Display','final','Algorithm','sqp' ,'MaxFunEvals', 10e5,'PlotFcns',@optimplotfval,'MaxIter',1000,'TolX',1e-10);
%options=optimset('Display','final','Algorithm','sqp','PlotFcns',@optimplotfval);

lb=zeros(size(control_s));
ub=ones(size(control_s));

for i=1:s;
   ub(i,:)= cs(i,1).*maxsept;
end

ub(s+1,:)=shellmax

const=@(x)septics(x,s,cs) %contrain abatement to septic totals?

[x,fval,exitflag] =fmincon(mnpv,control_s,[],[],[],[],lb,ub,const,options);

[npvop damageop INop netop costop kgsop kgiop scapop]=estuaryi5(x,cs,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma);

%plots control and in estuary N over time
scapop(s+1,taum+1:end)=x(s+1,:);

subplot(2,1,1)
plot(scapop(:,taum:end)')
subplot(2,1,2)
plot(INop(:,taum:end)-Target)

%plot(damage(1+taum:end))
