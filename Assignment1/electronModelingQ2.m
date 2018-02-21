clearvars
clearvars -Global
close all
format shorte

global C 
global Vx Vy x y xp yp xi yi Vxi Vyi Collisions
global numElect MarkerSize
global Mass T SavePics

numElect = 10000;
SavePics = 1;        %used at the end to save graphs on a 1, and not on a 0
numSteps = 1000;

len = 200e-9;
wid = 100e-9;

wallWidth = 1.2e-7;
wallX = .8e-7;
wallH1 = .4e-7;
wallH2 = .6e-7;

C.Mo = 9.10938215e-31;      % electron mass
C.kb = 1.3806504e-23;       % Blotzmann Const

T = 300;
Mass = 0.26*C.Mo;
k = 1.381 * 10 ^-23;
vth = sqrt(2*(C.kb*T)/(Mass));   %vth = 1.8702e5
dt = 10e-15;                     %10fs
TStop = numSteps*dt;

Prob = 1 - exp(-dt/.2e-12);   %probbility to interact with the backgorund

Limits = [0 len 0 wid];       %the drawing limits of the material simulated
MarkerSize = 1;

for i = 1:numElect            %initialize the position of each electron
                              %inside the material.
    x(i) = rand()*len;               
    y(i) = rand()*wid;
end

xi = x;
yi = y;

%averaging the distance to each neibouring
%electron to calculate mean free path
avgDistX(1:numElect) = x.*x(1:numElect);
avgDistY(1:numElect) = y.*y(1:numElect);    

Collisions = zeros(1,numElect);

xp = zeros(numElect);           %previous values will be used to track 
yp = zeros(numElect);           %the trajectories of the electrons

Vx = vth .* cos(2*pi*randn(1,numElect)); %initial velocities
Vy = vth .* sin(2*pi*randn(1,numElect));    
                                         
Vxi = Vx;
Vyi = Vy;

Vt = sqrt(Vx.*Vx + Vy.*Vy);
avgVel = sum(Vt)/numElect;

%histogram
figure(2);
histogram(Vt,200);
title('Average Thermalized Velocities');
xlabel('Thermal Velocity (m/s)');
ylabel('Amount per Bin');

fprintf('The Avg velocity is: %e; Calculated Thermal Velocity: %e \n'...
    , avgVel, vth);   

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
figure(1);
subplot(2,1,2);                     
axis([0 TStop 0 400]);
title('Material Temperature');
xlabel('Time (seconds)');
ylabel('Temp (Kelvin)');
hold on;
grid on;

%Find the initial temp of the material
for i = 1:numElect                  
    tempSum  = tempSum + (Mass*Vt(i)^2)/(2*C.kb);
end

avgTemp = tempSum/numElect;
Temp = [300 avgTemp];
Time = [0 t];
plot(t, avgTemp, '-');

numVisable = 10;                %This sets the amount of visable electrons
colorVec = hsv(numVisable);
tempSum = 0;                    %Reseting some values to zero to ensure
avgTemp = 0;                    %proper calculations
Vt = 0;
prevTemp = 0;

sumCollision = 0;               %initializing some helpers to calculate
sumCollTime = 0;                %the average collision time
numColl = 0;

while t < TStop                 %Loop to calcualte pos, and temp
    xp(1:numElect) = x(1:numElect);
    yp(1:numElect) = y(1:numElect);
    
    %update position before the bounds check
    %the bounds will rewrite this if an electron
    %is outside the bounds
    x(1:numElect) = x(1:numElect) + (dt .* Vx(1:numElect)); 
    y(1:numElect) = y(1:numElect) + (dt .* Vy(1:numElect)); 
    
    for i=1:numElect            %Loop to calcuate the boundaries, left and 
                                %right are periodic, the top and bottom
                                %are reflections
                                
       %Boundary conditions, not rethermalized
       if x(i) >= len
           xp(i) = 0;
           x(i) = dt * Vx(i);
       end
       if x(i) <= 0
           xp(i) =  xp(i) + len;
           x(i) =  xp(i) + dt*Vx(i);
       end   
       if y(i) >= wid || y(i) <= 0
           Vy(i) = - Vy(i);
       end
       
       %implement scattering here, the velocity is re-thermalized
       if rand() < Prob
          Vx(i) = vth * cos(2*pi*randn());
          Vy(i) = vth * sin(2*pi*randn());
          
          sumCollTime = sumCollTime + Collisions(i);%take the time of the 
          Collisions(i) = 0;             %last walk, reset the time between 
          numColl = numColl + 1;         %collisions, count the number of 
       end                               %collisions
       %sum the time between collions per electron
       Collisions(i) = Collisions(i) + dt;          
       
       Vt = sqrt(Vx(i)^2 + Vy(i)^2);            %As we loop to check bounds
       tempSum = tempSum + (Mass*Vt^2)/(2*C.kb);%we might aswell do the 
                                                %temp cacluations
       
       if i <= numVisable             %plot the difference in position,
           figure(1);                 %but only a small number will show
           subplot(2,1,1);
           plot([xp(i) x(i)], [yp(i) y(i)],'color',colorVec(i,:));
       end
    end
   
    avgTemp = tempSum/numElect;       %evaluate the avg temp of the system
    Temp = [prevTemp avgTemp];        %takes two points to make a line
    Time = [(t-dt) t];                %the previous temp and the previous
    figure(1);                        %time should line up, so t-dt is the
    subplot(2,1,2);                   %previous temp                                       
    plot(Time, Temp, '-', 'color', colorVec(1,:));
        
    prevTemp = avgTemp;               %used to calculate the material temp
    avgTemp = 0;
    tempSum = 0;
    pause(0.00001);
    t = t + dt;
    hold on;
end



%mean free path calculation
avgx = sum(Vxi - Vx);
avgy = sum(Vyi - Vy);
AvgDist = sum(sqrt(avgDistX.^2 + avgDistY.^2))/numElect;
avgTot = sqrt(avgx^2 + avgy^2)/sqrt(2)*pi*numElect*(AvgDist)^2;

%mean time between collisions
avgCollTime = sumCollTime / numColl;

fprintf('Mean Free Path Calcuated: %g Avg time between collisions: %g\n'...
    , avgTot, avgCollTime); 

%save the final state of the graphs to add to the report
if SavePics
    figure(1);
    saveas(gcf, 'ElectronsInSiliconQ2.jpg');

    figure(2);
    saveas(gcf, 'VelocityHistQ2.jpg');
end