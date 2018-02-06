clearvars
clearvars -GLOBAL
close all
global V;

N = 100;
M = 100;
R = 1;
L = 0;
T = 0;
B = 0;

V = rand(N,M);
numSims = 1000;

figure(1);
Ex = M;
Ey = N;

for i =1:numSims
    PA3(N,M,R,L,T,B,0,0); 
    figure(1);
    surf(1:M,1:N,V);
    fprintf("percent complete: %g \n", i/numSims *100);
    [Ex, Ey] = gradient(V);
    figure(2);
    surf(1:M,1:N,Ex);
    figure(3);
    surf(1:M,1:N,Ey);
    figure(4);
    quiver(Ex,Ey);
    pause(0.00000001);

end

