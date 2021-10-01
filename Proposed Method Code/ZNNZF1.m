function dot_u=ZNNZF1(t,u,k)
n=3;
m=2;
p=5;
h=10;
gamma=20;
delta=0.0001*ones(p,1);
x=u(1:n);
mu=u(n+m+2:n+m+p+1);
mB=[1 0 -1;0 1 -1;-1 0 -1;0 -1 -1;0 0 1];
b2=[0;0;0;0;h];
e=b2-mB*x;
D1=diag(e./sqrt(e.*e+mu.*mu+delta));
D2=diag(mu./sqrt(e.*e+mu.*mu+delta));

syms s1
A=diff(mA(s1));
s1=t;
dA=eval(A);
syms s2
b=diff(b1(s2));
s2=t;
db=eval(b);


C=[0;0;1];
mW=[0,0,0;0,0,0;0,0,0];

z1=zeros(n,n);
z2=zeros(n,p);
z3=zeros(p,n);
z4=zeros(p,p);
z5=zeros(n,1);
z6=zeros(p,1);
I=eye(p,p);
I2=eye(m+n+p+1,m+n+p+1);

mH=[mW mA(t)' mB';mA(t) z1 z2;-mB z3 I];
mM=[mW mA(t)' mB';mA(t) z1 z2;(D1-I)*mB z3 I-D2]+0.0001*I2;
mN=[z1 dA' z2;dA z1 z2;z3 z3 z4];

qt=[C;-b1(t);b2-sqrt(e.*e+mu.*mu+delta)];
dqt=[z5;-db;z6];
m_M=inv(mM);
err=mH*u(1:m+1+n+p)+qt;

if k==2
    derr=-gamma*err;
    dot_u=m_M*(-mN*u(1:m+1+n+p)-dqt+derr); 
end 

if k==3
     dot_u=m_M*(-mN*u(1:m+1+n+p)-dqt-gamma*SSBP(mH*u(1:m+1+n+p)+qt));
end
if k==4
    derr=-gamma*sat(err,gamma);
    dot_u=m_M*(-mN*u(1:m+1+n+p)-dqt+derr); 
end 
result=dot_u; 

t