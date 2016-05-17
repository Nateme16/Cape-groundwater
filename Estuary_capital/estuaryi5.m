function [npv damage IN net cost kgs kgi scap] = estuaryi5(control_s,cs,taum,n,k,tau,s,r,cost_s,cost_in,R,IN0,B,Target,gamma)
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

for ii=1:s; %pre-control period attenuated controllable load
   cstart(ii)=cs(ii,1)*exp(-k*tau(ii));
end

TCstart=sum(sum(cstart)); 
TC=zeros(1,n);
TC(:,:)=TCstart;

IN=zeros(1,n);
IN(:,1)=IN0;

for i=2:taum; %creates pre-control N stock in estuary, "do nothing scenario"
    IN(:,i)=dettman(TC(i),R,IN(:,i-1),B); 
end

scap=zeros(s,n);

for i=taum+1:taum+n; % this loop runs the control through time
 
    for ii=1:s; 
    scap(ii,i)=scap(ii,i-1) + control_s1(ii,i); %source tech capital
    c2(ii,1)=(cs(ii,i-tau(ii)) - scap(ii,i-tau(ii)) )*exp(-k*tau(ii));%calculates total load depending on year and control with delay time tau and attenuation
    end
    
    tc1=sum(sum(c2)); % sums total load
    
    IN1(i)=dettman(tc1,R,IN(i-1),B); % calculates N stock for year i

    IN(i)= IN1(i)-control_s1(s+1,i); % calculates N stock for year i after in estuary treatment
    
    
cost(i)=sum(cost_s.*(scap(:,i))) + (cost_in.*(control_s1(s+1,i).^2)); %creates cost of year i based on choices in control matrix

kgs(i)=sum(scap(:,i));
kgi(i)=control_s1(s+1,i);

damage(i)=((IN(i)-Target)^2)*gamma; %damages of N stock in estuary

net(i)=(damage(i)+cost(i))*exp(-r*i); %  damages+costs discounted

end
npv=sum(net);

end