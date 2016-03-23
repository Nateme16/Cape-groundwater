function [npv damage IN net] = estuaryi2(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B)
% esturaryi2  - Simulates loading and estuary N stock
% control_s - control matrix
% c - source loadings uncontrolled
% taum - the max of the delay times (this is to deal with the pre-control
% time
% n - number of years
% k - % N removed per year in transit
% tau - matrix of time delay for sources
% s- number of sources
% cost_s- cost parameter for sources
% cost_in - cost parameter for in-estuary treatment
% R - net N caryover from year before
% IN0 - starting N stock
% B - background annual loading from non-controllable sources

c2=zeros(1,s)';
cost=zeros(1,s)';
control_s1=[ones(s+1,taum) control_s];
TC=sum(c,1); 
IN=zeros(1,n);
IN(:,1)=IN0;

for i=2:taum; %creates pre-control N stock
    IN(:,i)=dettman(TC(i),R,IN(:,i-1),B); 
end

for i=taum+1:taum+n; % this loop runs the control through time
 
    for ii=1:s; %calculates total load depending on year with delay time tau
    c2(ii,1)=(control_s1(ii,i-tau(ii))*c(ii,i-tau(ii)))*exp(-k*tau(ii));
    end
    
    tc1=sum(sum(c2)); % sums total load
    
    IN(i)=dettman(tc1,R,TC(i-1),B); % calculates N stock for year i

    TC(i)= IN(i).*control_s1(s+1,i); % calculates N stock for year i after in estuary treatment
    
    
cost(i)=sum(cost_s.*(1./control_s1(1:s,i)-1).*c(1:s,i)) + (cost_in.*(1./control_s1(s+1,i)-1)).*IN(i); %creates cost of year i based on choices in control matrix

damage(i)=(100*TC(i)^2); %damages of N stock in estuary

net(i)=(damage(i)+cost(i))*exp(-r*i); %  damages+costs

end

%This part was when I was using fminsearch, needed to penalize non-feasible
%options
if control_s1(:,:)>0 & control_s1(:,:)<=1;
npv=sum(net);

else
    npv=10^100;
end
end