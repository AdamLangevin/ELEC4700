clearvars
clearvars -Global
close all
format shorte

global C 
global Vx Vy x y xp yp
global numElect MarkerSize
global Mass T SavePics

numElect = 10000;
SavePics = 1;
numVisable = 10;         %This sets the amount of visable electrons
numSteps = 1000;

len = 200e-9;
wid = 100e-9;

C.Mo = 9.10938215e-31;      % electron mass
C.kb = 1.3806504e-23;       % Blotzmann Const

T = 300;
Mass = 0.26*C.Mo;
k = 1.381 * 10 ^-23;
vth = sqrt(2*(C.kb*T)/(Mass));
dt = 10e-15;
TStop = numSteps*dt;
Limits = [0 len 0 wid];
MarkerSize = 1;

%initialize the position of each electron
%inside the material.
for i = 1:numElect              
    x(i) = rand()*len;          
    y(i) = rand()*wid;
end

%previous values will be used to track 
%the trajectories of the electrons

xp = zeros(numElect);               
yp = zeros(numElect); 

%initial velocities
Vx(1:numElect) = vth * cos(2*pi*randn(1,numElect)); 
Vy(1:numElect) = vth * sin(2*pi*randn(1,numElect)); 

Vt = sqrt(Vx.*Vx + Vy.*Vy);
tempSum = 0;

t = 0;

%initialize the electron position plot
figure(1);                       
subplot(2,1,1);
axis(Limits);
title('Electron Movement Through Silicon');
xlabel('X');
ylabel('Y');
hold on;
grid on;

%initialize the material temperuature plot
subplot(2,1,2);                  
axis([0 TStop 0 400]);
title('Material Temperature');
xlabel('Time (seconds)');
ylabel('Temp (Kelvin)');
hold on;
grid on;
        
for i = 1:numElect         %Find the initial temp of the material
    tempSum  = tempSum + (Mass*Vt(i)^2)/(2*C.kb);
end
avgTemp = tempSum/numElect;
Temp = [300 avgTemp];
Time = [0 t];
plot(t, avgTemp, '-');

colorVec = hsv(numVisable);      %Random color assignments
tempSum = 0;                     %Reseting some values to zero 
avgTemp = 0;                     %to ensure proper calculations
Vt = 0;
prevTemp = 0;

while t < TStop                  %Loop to calcualte pos, and temp
    xp = x;
    yp = y;
    
    x(1:numElect) = x(1:numElect) + (dt .* Vx(1:numElect));
    y(1:numElect) = y(1:numElect) + (dt .* Vy(1:numElect));
    
    for i=1:numElect  %Loop to calcuate the boundaries, left and 
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
              
       Vt = sqrt(Vx(i)^2 + Vy(i)^2);       %As we loop to check bounds
       tempSum = tempSum + (Mass*Vt^2)/(2*C.kb);%we might aswell do 
                                           %the temp cacluations
       
       X = [xp(i) x(i)];
       Y = [yp(i) y(i)];   %reduce this to the inside of the plot
       if i < numVisable
           subplot(2,1,1);
           plot(X,Y,'color',colorVec(i,:));
       end
    end
   
    avgTemp = tempSum/numElect;%evaluate the avg temp of the system
    Temp = [prevTemp avgTemp]; %takes two points to make a line
    Time = [(t-dt) t];         %the previous temp and the previous
    subplot(2,1,2);            %time should line up, so t-dt is the
                               %previous temp
    plot(Time, Temp, '-', 'color', colorVec(1,:));
        
    prevTemp = avgTemp;                 
    avgTemp = 0;
    tempSum = 0;
    pause(0.00001);
    t = t + dt;
end

if SavePics
    figure(1);
    saveas(gcf, 'ElectronsInSiliconQ1.jpg');
end
