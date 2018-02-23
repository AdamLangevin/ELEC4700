clear
clf

n = 100;   %no of points
dx = 1/n;  %step size
G = zeros(n-2);  %initialize G matrix
V = zeros(n+1,1); %initial solution of potential

for i=1:n-1
   G(i,i) = -2; %cofficient for 1D finite difference eqn
end

for i=1:n-2
   G(i,i+1) = 1; %cofficient for 1D finite difference eqn
   G(i+1,i) = 1; %cofficient for 1D finite difference eqn
end

B = zeros(n-1,1); %initialize B matrix
v0 = 1; %right BC
v1 = 0; %left BC

B(1) = v0;   %right BC
B(n-1) = v1; %left BC

V = (1/(dx.^2)*G)\(-B*(1/dx.^2)); %solution of potential
x = linspace(0,1,n-1); %range to plot
figure(6);
plot(x,V)
title('Potential vs. Position')
xlabel('X Position')
ylabel('Potential (V)') 