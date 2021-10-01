function output=sat(E,gamma,r)
if nargin==2,r=0.8;
end
for i=1:length(E)    
    temp=sign(E(i))*gamma*(abs(E(i)).^r);   
    if temp<-1
        E(i)=-1;
    else
        if temp<=1
           E(i)=temp;
        else
           E(i)=1;
        end
    end
end
output=E;


