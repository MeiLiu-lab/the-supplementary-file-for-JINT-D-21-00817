clc;
clear;
close all;

m=2;n=3;p=5;
T=10;
u0=0.0001*rand(m+n+p+1,1);
tspan=linspace(0,T,600);
options=odeset(); 

[t,u1(:,:)]=ode45(@ZNNZF1rightside,tspan,u0,options,1);%AFznn
[t,u2(:,:)]=ode45(@ZNNZF1,tspan,u0,options,2);%znn
[t,u3(:,:)]=ode45(@ZNNZF1,tspan,u0,options,3);%SSBP
[t,u4(:,:)]=ode45(@ZNNZF1,tspan,u0,options,4);%satZNN

Num1=length(t);
for i=1:Num1
    [mH1,qt1]=solveHv(t(i),u1(i,1:m+n+p+1));
     NORM1(i)=norm(mH1*u1(i,1:m+n+p+1)'+qt1);
    
    [mH2,qt2]=qiujie(t(i),u2(i,1:m+n+p+1));
    NORM2(i)=norm(mH2*u2(i,1:m+n+p+1)'+qt2); 
    
    [mH3,qt3]=qiujie(t(i),u3(i,1:m+n+p+1));
    NORM3(i)=norm(mH3*u3(i,1:m+n+p+1)'+qt3); 
    
    [mH4,qt4]=qiujie(t(i),u4(i,1:m+n+p+1));
    NORM4(i)=norm(mH4*u4(i,1:m+n+p+1)'+qt4);  
    
    
end
figure (1);
plot(t,NORM1(:),t,NORM2(:),t,NORM3(:),t,NORM4(:),'MarkerSize',2.2,'linewidth',1.1);

