clear all
clc
%Diffusion Equation
%Explicit
%%Bc1-1
N=15;
k=1;
%%Distance mesh y and x
y=linspace(-pi,pi,N);            %Evenly space y over N increments
x=linspace(-pi,pi,N);            %Evenly space x over N increments
dy=y(N)-y(N-1);
dx=x(N)-x(N-1);
dt=dx^2/(6*k);                   %Stable is dt<0.5*dx^2/k

%%Define constants and given equations
ax=-pi;                                  %Left side of square 
bx=pi;                                   %Right side of square
ay=-pi;                                  %Bottom of square
by=pi;                                   %Top of square
fb=(bx-x).^2.*cos(pi*x/bx);              %Top BC
gb=x.*(bx-x).^2;                         %Bottom BC
fbax=(bx-ax)^2*cos(pi*ax/bx);            %gb(ax)
gbax=ax*(bx-ax)^2;                       %fb(ax)
u_ax=(gbax+(y-ay)./(by-ay)*(fbax-gbax)); %Left BC
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
    
for j=2:N-1
   for i=2:N-1
    u_t(i,j)=r*u(i,j)+lam*u(i+1,j)+lam*u(i,j+1)+lam*u(i-1,j)+lam*u(i,j-1);
   end
end
%%Boundary Conditions
u_t(1,:)=fb;                            %Apply top BC to u at time+dt
u_t(N,:)=gb;                            %Apply bottom BC to u at time+dt
u_t(:,1)=u_ax;                          %Apply left BC to u at time+dt

for i=2:N-1
    u_t(i,N)=r*u(i,j)+lam*u(i,j+1)+2*lam*u(i-1,j)+lam*u(i,j-1);
end
    u=u_t;
    %Calculate error
    err_1=sum(abs(u_t-u_1));
    err=sum(err_1);
    
  
u_1=u_t;                                %Old version of u=new version of u
iteration=iteration+1;
t=t+dt;
%Create countour plot of u throughout time until steady state is reached
drawnow                                 
contourf(X,Y,u)
caxis([min(min(u)),max(max(u))])
colorbar
title(sprintf('%11.3f  seconds',t))
xlabel('x axis')
ylabel('y axis')

end
sprintf('Iterations %d',iteration)

    
