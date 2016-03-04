function [npv damage TC net ] = estuarystatic(control_s,c,taum,n,k,tau,s,r,cost_s)

c2=zeros(1,s)';
cost=zeros(1,s)';
control_s1=[ones(s,taum) repmat(control_s,1,n)];
TC= sum(c.*repmat(exp(-k*tau),1,size(c,2),1));

for i=taum+1:taum+n;
 
    for ii=1:s;
    c2(ii,1)=(control_s1(ii,i-tau(ii))*c(ii,i-tau(ii)))*exp(-k*tau(ii));
    end
    
    TC(i)= sum(sum( c2 ));
    
cost(i)=sum(cost_s.*((1./control_s1(:,i)-1).*c(:,i)));

damage(i)=(100*TC(i)^2);

net(i)=(damage(i)+cost(i))*exp(-r*i);

end


if control_s1(:,:)>0 & control_s1(:,:)<=1
npv=sum(net);

else
    npv=10^100;
end
end
