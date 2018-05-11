clear all
clc
%Diffusion Equation
%ADI
%%Bc1-1
%%Discretize x and y axes and define space and time steps
space=15;
n=30;
k=1;
y=linspace(-pi,pi,space);  %
x=linspace(-pi,pi,space);
dy=y(space)-y(space-1);
dx=x(space)-x(space-1);
dt=dx^2/(6*k);
[X,Y]=meshgrid(x,y);
%%Define given equations and constants
ax=-pi;                                   %Left side of square
bx=pi;                                    %Right side of square
ay=-pi;                                   %Bottom of square
by=pi;                                    %Top of square
fb=(bx-x).^2.*cos(pi*x/bx);               %Top BC
gb=x.*(bx-x).^2;                          %Bottom BC
fbax=(bx-ax)^2*cos(pi*ax/bx);             %fb(ax)
gbax=ax*(bx-ax)^2;                        %gb(ax)
u_ax=(gbax+(y-ay)./(by-ay)*(fbax-gbax));  %Left BC
Neumann=zeros(1,length(x));               %Right BC

N=length(x)-2;                            %Define length used ADI scheme

lam=dt/dx^2;
L=2*(1-lam);                                

%Define tri-diagonal matrices
m1=2*(1+lam)*ones(1,N); 
m2=-lam*ones(1,N-1);
M=diag(m1);
M=diag(m2,1)+M;
M=diag(m2,-1)+M;                          %Used to solve for u
m=M;                                      
m(N,N+1)=M(1,2);
m(N+1,[N,N+1])=[-2*lam,M(N,N)];           %Used to solve for u and Neumann BCs

%%Solve for u
%Assign Boundary Conditions
u1=zeros(space,space);
u1(1,:) = fb;               %fb=Top BC
u1(space,:) =gb;            %gb=Bottom Bc, 
u1(1:1:N+2,1)=u_ax;        %Left BC
u0=u1;
u=zeros(space-2,space-2);

iteration=0;
t=0;
err=10000;
while err>0.0001
    
    %Sweep in y direction
       %Column 1
       B(1:N,1)=lam*u_ax(2:N+1)'+L*u(1:N,1)+lam*u(1:N,2);
       B(1)=B(1)+lam*gb(2);
       B(N)=B(N)+lam*fb(2);
       u(:,1)=M\B;
       
       %Columns 2 through N-1
       B2(1:N,1:N-2)=lam*u(1:N,1:N-2)+L*u(1:N,2:N-1)+lam*(u(1:N,3:N));
       B2(1,:)=B2(1,:)+lam*gb(3:N);
       B2(N,:)=B2(N,:)+lam*fb(3:N);
       u(:,2:N-1)=M\B2;
             
       %Column N
       B(1:N,1)=lam*u(1:N,N-1)+L*u(1:N,N)+lam*Neumann(2:N+1)';
       B(1)=B(1)+lam*gb(N+1);
       B(N)=B(N)+lam*fb(N+1);
       u(:,N)=M\B;
       
       %Add Neumann BCs to x=ay
       B(1:N,1)=(2*lam)*u(1:N,N)+L*Neumann(2:N+1)';
       B(1)=B(1)+lam*gb(N+1);
       B(N)=B(N)+lam*fb(N+1);
       Neumann(2:N+1)=(M\B)';
      
      %Sweep in x direction
       %Row 1
       b(1:N,1)=lam*gb(2:N+1)+L*u(1,1:N)+lam*u(2,1:N);
       b(1)=b(1)+lam*u_ax(2);
       b(N+1)=lam*gb(N+2)+L*Neumann(2)+lam*Neumann(3);
       u2=m\b;
       u(1,:)=u2(1:N);
       Neumann(2)=u2(N+1);
       
       %Rows 2 through N-1
       b2(1:N-2,1:N)=lam*u(1:N-2,1:N)+L*u(2:N-1,1:N)+lam*u(3:N,1:N);
       b2(:,1)=b2(:,1)+lam*u_ax(3:N)';
       b2(:,N+1)=lam*Neumann(2:N-1)+L*Neumann(3:N)+lam*Neumann(4:N+1);
       u2=(m\b2')';
       u(2:N-1,:)=u2(1:N-2,1:N);
       Neumann(3:N)=u2(:,N+1);
       
       %Row N
       b(1:N,1)=lam*u(N-1,1:N)+L*u(N,1:N)+lam*fb(2:N+1);
       b(1)=b(1)+lam*u_ax(N+1);
       b(N+1)=lam*Neumann(N)+L*Neumann(N+1)+lam*fb(N+2);
       u2=m\b;
       u(N,:)=u2(1:N);
       Neumann(N+1)=u2(N+1);
       
       u1(1:1:N+2,space)=Neumann;
       u1(2:1:N+1,2:N+1)=u;
       
       
       
       err_1=sum(abs(u1-u0));
       err=sum(err_1);
       
       u0=u1;
       iteration=iteration+1;
       %Create a countour plot of u throughout time until steady state is
       %reached
       drawnow
       t=t+dt;
       drawnow
       contourf(X,Y,u1)
       caxis([min(min(u1)),max(max(u1))])
        colorbar
       title(sprintf('%11.3f  seconds',t))
       xlabel('x axis')
       ylabel('y axis')
       
end


sprintf('Iterations %d',iteration)

