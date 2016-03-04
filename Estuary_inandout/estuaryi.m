function [npv damage TC net] = estuaryi(control_s,c,taum,n,k,tau,s,r,cost_s,cost_in)

c2=zeros(1,s)';
cost=zeros(1,s)';
control_s1=[ones(s+1,taum) control_s];
TC=sum(c,1);

for i=taum+1:taum+n;
 
    for ii=1:s;
    c2(ii,1)=(control_s1(ii,i-tau(ii))*c(ii,i-tau(ii)))*exp(-k*tau(ii));
    end
    
    TC1(i)=sum(sum( c2 ));
    TC(i)= sum(sum( c2 )).*control_s1(s+1,i);
    
cost(i)=sum(cost_s.*(1./control_s1(1:s,i)-1).*c(1:s,i)) + (cost_in.*(1./control_s1(s+1,i)-1)).*TC1(i);

damage(i)=(100*TC(i)^2);

net(i)=(damage(i)+cost(i))*exp(-r*i);

end


if control_s1(:,:)>0 & control_s1(:,:)<=1;
npv=sum(net);

else
    npv=10^100;
end
end