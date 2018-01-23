clearvars
clearvars -Global
close all
format shorte

global C 
global Vx Vy x y
global numElect MarkerSize
global Mass T

numElect = 1000;

len = 200e-9;
wid = 100e-9;

C.Mo = 9.10938215e-31;      % electron mass
C.kb = 1.3806504e-23;       % Blotzmann Const

T = 30;
Mass = 0.26*C.Mo;
k = 1.381 * 10 ^-23;
vth = sqrt(2*(C.kb*T)/(Mass));
dt = 10e-15;
TStop = 1000*dt;

Limits = [0 len 0 wid];
MarkerSize = 1;

for i = 1:numElect
    x(i) = rand()*200e-9;
    y(i) = rand()*100e-9;
end

Vx(1:numElect) = vth * cos(2*pi*randn(1,numElect));
Vy(1:numElect) = vth * sin(2*pi*randn(1,numElect));

t = 0;
while t < TStop
    x(1:numElect) = x(1:numElect) + (dt .* Vx(1:numElect));
    y(1:numElect) = y(1:numElect) + (dt .* Vy(1:numElect));
    
    t = t + dt;
    for i=1:numElect
       if x(i) >= len
           x(i) = x(i) - len;
       end
       if x(i) <= 0
           x(i) = x(i) + len;
       end   
       if y(i) >= wid || y(i) <= 0
           Vy(i) = - Vy(i);
       end
    end

    fprintf('time: %g (%5.2g %%)\n', t, t/TStop*100);
    numVisable = 20;
    for j=1:numVisable
        colorVec = hsv(numVisable +1);
        plot(x(j),y(j), 'o', 'markers', MarkerSize, 'color', colorVec(j,:), 'MarkerFaceColor', colorVec(j,:));
    end
    axis(Limits);
    hold on;
    pause(0.00001);
end


