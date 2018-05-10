%%Distance mesh y and x
space=40;
n=20;
k=1;
y=linspace(-pi,pi,space);
x=linspace(-pi,pi,space);
dy=y(space)-y(space-1);
dx=x(space)-x(space-1);
dt=dx^2/(6*k);
%%Define Initial and Boundary Conditions
ax=-pi;
bx=pi;
ay=-pi;
by=pi;
fb=(bx-x).^2.*cos(pi*x/bx);               %Bottom BC
gb=x.*(bx-x).^2;                          %Top BC
fbax=(bx-ax)^2*cos(pi*ax/bx);
gbax=ax*(bx-ax)^2;
u_ax=(gbax+(y-ay)./(by-ay)*(fbax-gbax));  %Left BC
Neumann=zeros(1,length(x));              %Right BC

N=length(x)-2;
Nt=2:N+1;
lam=dt/dx^2;
L=(1-lam)*2;

m1=ones(1,N)*2*(1+lam);
m2=-lam*ones(1,N-1);
M=diag(m1);
M=diag(m2,1)+M;
M=diag(m2,-1)+M;
m=M;
m(N,N+1)=M(1,2);
m(N+1,[N,N+1])=[-2*lam,M(N,N)];

%%Start solving for u
%Assign Boundary Conditions
u1=zeros(space,space);

u1(1:1:N+2,1)=u_ax;      %Left BC
u1([1,space],:) = [gb;fb];
u=zeros(space-2,space-2);
iteration=0;

while iteration<10^4
    
    %Sweep in y direction
       %Column 1
       B(1:N,1)=lam*u_ax(Nt)'+L*u(1:N,1)+lam*u(1:N,2);
       B(1)=B(1)+lam*fb(2);
       B(N)=B(N)+lam*gb(2);
       u(:,1)=M\B;
       
       %Columns 2 through N-1
       B2(1:N,1:N-2)=L*u(1:N,2:N-1)+lam*(u(1:N,3:N))+lam*u(1:N,1:N-2);
       B2(1,:)=B2(1,:)+lam*fb(3:N);
       B2(N,:)=B2(N,:)+lam*gb(3:N);
       u(:,2:N-1)=M\B2;
             
       %Column N
       B(1:N,1)=L*u(1:N,N)+lam*u(1:N,N-1)+lam*Neumann(2:N+1);
       B(1)=B(1)+lam*u_ax(N+1);
       B(N)=B(N)+lam*gb(N+1);
       u(:,N)=A\B';
       
        iteration=iteration+1;
end




