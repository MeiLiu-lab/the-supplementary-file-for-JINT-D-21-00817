function F=F1(x,t)
    tc=0.5;
for i=1:length(x)
     F =(exp(i)-1)/((tc-t)*(exp(i)));  
end