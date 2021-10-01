function [mH,qt]=qiujie(t,u)
n=3;
m=2;
p=5;
h=10;
delta=0.0001*ones(p,1);

x=u(1:n)';
mu=u(n+m+2:n+m+p+1)';
mB=[1 0 -1;
    0 1 -1;
    -1 0 -1;
    0 -1 -1
    0 0 1];
b2=[0;0;0;0;h];
e=b2-mB*x;
C=[0;0;1];
% b1=[sin(4*t);sin(4*t);sin(4*t)];
% mA=[sin(4*t),cos(4*t),0;
%     sin(4*t),cos(4*t),0;
%     sin(4*t),cos(4*t),0];
mW=[0,0,0;0,0,0;0,0,0];
% z1=zeros(n,n);
% z2=zeros(m,p);
% z3=zeros(p,m);
% z4=zeros(p,p);
% z5=zeros(m,1);
z1=zeros(n,n);
z2=zeros(n,p);
z3=zeros(p,n);
% z4=zeros(p,p);
% z5=zeros(n,1);
I=eye(p,p);
mH=[mW mA(t)' mB';mA(t) z1 z2;-mB z3 I];
qt=[C;-b1(t);b2-sqrt(e.*e+mu.*mu+delta)];
t


