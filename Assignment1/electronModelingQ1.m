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
%Lambda = log(1+ Prob)              %mean free path?
Limits = [0 len 0 wid];
MarkerSize = 1;

for i = 1:numElect                  %initialize the position of each electron
    x(i) = rand()*200e-9;           %inside the material. 
    y(i) = rand()*100e-9;
end

xp = zeros(numElect);               %previous values will be used to track 
yp = zeros(numElect);               %the trajectories of the electrons

Vx(1:numElect) = vth * cos(2*pi*randn(1,numElect)); %initial velocities
Vy(1:numElect) = vth * sin(2*pi*randn(1,numElect)); %based on the thermal velocity

Vt = sqrt(Vx.*Vx + Vy.*Vy);
tempSum = 0;

t = 0;

figure(1);                          %initialize the electron position plot
subplot(2,1,1);
axis(Limits);
title('Electron Movement Through Silicon');
xlabel('X');
ylabel('Y');
hold on;
grid on;

subplot(2,1,2);                     %initialize the material temperuature plot
axis([0 TStop 0 400]);
title('Material Temperature');
xlabel('Time (seconds)');
ylabel('Temp (Kelvin)');
hold on;
grid on;
        
for i = 1:numElect                  %Find the initial temp of the material
    tempSum  = tempSum + (Mass*Vt(i)^2)/(2*C.kb);
end
avgTemp = tempSum/numElect;
Temp = [300 avgTemp];
Time = [0 t];
plot(t, avgTemp, '-');

numVisable = 10;                    %This sets the amount of visable electrons
colorVec = hsv(numVisable + 1);     %and adds different color values to each vector

tempSum = 0;                        %Reseting some values to zero to ensure
avgTemp = 0;                        %proper calculations
Vt = 0;
prevTemp = 0;

while t < TStop                     %Loop to calcualte pos, and temp
    xp(1:numElect) = x(1:numElect);
    yp(1:numElect) = y(1:numElect);
    
    x(1:numElect) = x(1:numElect) + (dt .* Vx(1:numElect));
    y(1:numElect) = y(1:numElect) + (dt .* Vy(1:numElect));
    
    for i=1:numElect                %Loop to calcuate the boundaries, left and 
                                    %right are periodic, the top and bottom
                                    %are reflections
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
       
       %implement scattering here
       
       Vt = sqrt(Vx(i)^2 + Vy(i)^2);                %As we loop to check bounds
       tempSum = tempSum + (Mass*Vt^2)/(2*C.kb);    %we might aswell do the temp
                                                    %cacluations
       
       X = [xp(i) x(i)];
       Y = [yp(i) y(i)];            %reduce this to the inside of the plot
       if i < numVisable
           subplot(2,1,1);
           plot(X,Y,'color',colorVec(i,:));
       end
    end
   
    avgTemp = tempSum/numElect;         %evaluate the avg temp of the system
    Temp = [prevTemp avgTemp];          %takes two points to make a line
    Time = [(t-dt) t];                  %the previous temp and the previous
    subplot(2,1,2);                     %time should line up, so t-dt is the
                                        %previous temp
    plot(Time, Temp, '-', 'color', colorVec(1,:));
    
    %fprintf('time: %g (%5.2g %%)\n', t, t/TStop*100);
    
    prevTemp = avgTemp;                 
    avgTemp = 0;
    tempSum = 0;
    pause(0.00001);
    t = t + dt;
end


