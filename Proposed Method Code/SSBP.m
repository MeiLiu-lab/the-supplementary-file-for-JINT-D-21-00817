function output=SSBP(X,r)
if nargin==1,r=0.8;
end
output=zeros(size(X));
for i=1:length(X)    
      output(i) =(abs(X(i)).^r).*sign(X(i));   
end
%output=0.5*(abs(X).^r).*sign(X)+0.5*(abs(X).^(1/r)).*sign(X);
