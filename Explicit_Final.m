clear all
clc



N=20;
k=1;
%%Distance mesh y and x
y=linspace(-pi,pi,N);
x=linspace(-pi,pi,N);
dy=y(N)-y(N-1);
dx=x(N)-x(N-1);
dt=dx^2/(6*k);                   %Stable is dt<0.5*dx^2/k

%%Initial and Boundary Conditions
ax=-pi;
bx=pi;
ay=-pi;
by=pi;
fb=(bx-x).^2.*cos(pi*x/bx);              %Bottom BC
gb=x.*(bx-x).^2;                         %Top BC
fbax=(bx-ax)^2*cos(pi*ax/bx);
gbax=ax*(bx-ax)^2;
u_ax=fliplr(gbax+(y-ay)./(by-ay)*(fbax-gbax));   %Left BC
%%Nested For Loops for Explicit Method
u_t=zeros(N,N);
lam=k*dt/(dx^2);
r=1-4*lam;
[X,Y]=meshgrid(x,y);
u=zeros(N,N);
t=0;                                    %Initial Time=0
err=10000;
iteration=0;
u_1=u;
while err>0.0001
    u_t=zeros(N,N);
for j=2:N-1
   for i=2:N-1
    u_t(i,j)=r*u(i,j)+lam*u(i+1,j)+lam*u(i,j+1)+lam*u(i-1,j)+lam*u(i,j-1);
   end
end

u_t(1,:)=gb;                            %Apply top BC to u at time+dt
u_t(N,:)=fb;                            %Apply bottom BC to u at time+dt
u_t(:,1)=u_ax;                          %Apply left BC to u at time+dt

for i=2:N-1
    u_t(i,N)=r*u(i,j)+lam*u(i,j+1)+2*lam*u(i-1,j)+lam*u(i,j-1);
end
    u=u_t;
    err_1=sum(abs(u_t-u_1));
    err=sum(err_1);
 u_1=u_t;
iteration=iteration+1;
t=t+dt;
drawnow
contourf(X,Y,u)
title(sprintf('%11.3f  seconds',t))
xlabel('x axis')
ylabel('y axis')
end
sprintf('Iterations %d',iteration)

    
