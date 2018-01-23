clearvars
clearvars -Global
close all
format shorte

global C 
global Vx Vy x y xp yp
global numElect MarkerSize
global Mass T

numElect = 1000;

len = 200e-9;
wid = 100e-9;

C.Mo = 9.10938215e-31;      % electron mass
C.kb = 1.3806504e-23;       % Blotzmann Const

T = 300;
Mass = 0.26*C.Mo;
k = 1.381 * 10 ^-23;
vth = sqrt(2*(C.kb*T)/(Mass));
dt = 10e-15;
TStop = 1000*dt;

%Prob = 1 - exp(-10e-15/.2e-12);    %probbility to interact with the backgorund
%log(1+ Prob)                       %mean free path?
Limits = [0 len 0 wid];
MarkerSize = 1;

for i = 1:numElect
    x(i) = rand()*200e-9;
    y(i) = rand()*100e-9;
end

xp = zeros(numElect);
yp = zeros(numElect);

Vx(1:numElect) = vth * cos(2*pi*randn(1,numElect));
Vy(1:numElect) = vth * sin(2*pi*randn(1,numElect));

figure(1);
axis(Limits);
hold on;
grid on;

numVisable = 20;
colorVec = hsv(numElect + 1);

t = 0;
while t < TStop
    xp(1:numElect) = x(1:numElect);
    yp(1:numElect) = y(1:numElect);
    
    x(1:numElect) = x(1:numElect) + (dt .* Vx(1:numElect));
    y(1:numElect) = y(1:numElect) + (dt .* Vy(1:numElect));
    
    for i=1:numElect
       if x(i) >= len
           xp(i) = 0;
           x(i) = dt * Vx(i);
       end
       if x(i) <= 0
           xp(i) = xp(i) + len;
           x(i) = xp(i) + dt*Vx(i);
       end   
       if y(i) >= wid || y(i) <= 0
           Vy(i) = - Vy(i);
       end
       
       %implment probability here?
       
       X = [xp(i) x(i)];
       Y = [yp(i) y(i)];
       plot(X,Y,'color',colorVec(i,:));
    end

    %fprintf('time: %g (%5.2g %%)\n', t, t/TStop*100);
    
    pause(0.00001);
    
    t = t + dt;
end


