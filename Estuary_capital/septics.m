function [c,ceq]=septics(x,s,cs)

for i=1:s
c(i)=  sum(x(i,:)) - cs(i,1) ;
end

ceq=[];
end
