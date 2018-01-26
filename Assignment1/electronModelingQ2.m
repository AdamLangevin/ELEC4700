clearvars
clearvars -Global
close all
format shorte

global C 
global Vx Vy x y xp yp xi yi Vxi Vyi Collisions
global numElect MarkerSize
global Mass T SavePics

numElect = 100;
SavePics = 1;

len = 200e-9;
wid = 100e-9;

C.Mo = 9.10938215e-31;      % electron mass
C.kb = 1.3806504e-23;       % Blotzmann Const

T = 300;
Mass = 0.26*C.Mo;
k = 1.381 * 10 ^-23;
vth = sqrt(2*(C.kb*T)/(Mass));      %vth = 1.8702e5
dt = 10e-15;                        %10fs
TStop = 1000*dt;

Prob = 1 - exp(-dt/.2e-12);         %probbility to interact with the backgorund
%L = log(1 - Prob);                 %mean free path?
%Lambda = L/diffL

Limits = [0 len 0 wid];
MarkerSize = 1;

for i = 1:numElect                  %initialize the position of each electron
    x(i) = rand()*len;           %inside the material. 
    y(i) = rand()*wid;
end

xi = x;
yi = y;

avgDistX(1:numElect) = x.*x(1:numElect);
avgDistY(1:numElect) = y.*y(1:numElect);

Collisions = zeros(1,numElect);

xp = zeros(numElect);               %previous values will be used to track 
yp = zeros(numElect);               %the trajectories of the electrons

Vx = vth .* cos(2*pi*randn(1,numElect)); %initial velocities
Vy = vth .* sin(2*pi*randn(1,numElect)); %based on the thermal velocity
                                         %should this be
                                         %sqrt(-log(numElect*sqrt(2*pi*C.kb*T/Mass))*2*C.kb*T/Mass)? 
                                         
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

%avgVel = 1.784874e+05
fprintf('The Average velocity is: %e; Calculated Thermal Velocity: %e \n', avgVel, vth);   

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

figure(1);
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
%colorVec = rand(numVisable,3);      %and adds different color values to each vector
colorVec = hsv(numVisable);
tempSum = 0;                        %Reseting some values to zero to ensure
avgTemp = 0;                        %proper calculations
Vt = 0;
prevTemp = 0;

sumCollision = 0;
sumCollTime = 0;
numColl = 0;

while t < TStop                     %Loop to calcualte pos, and temp
    xp(1:numElect) = x(1:numElect);
    yp(1:numElect) = y(1:numElect);
    
    x(1:numElect) = x(1:numElect) + (dt .* Vx(1:numElect));
    y(1:numElect) = y(1:numElect) + (dt .* Vy(1:numElect));
    
    for i=1:numElect                %Loop to calcuate the boundaries, left and 
                                    %right are periodic, the top and bottom
                                    %are reflections
       %Boundary hit conditions
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
          
          sumCollTime = sumCollTime + Collisions(i);%take the time of the last walk
          Collisions(i) = 0;                        %reset the time between collisions
          numColl = numColl + 1;                    %count the number of collisions
       end
       Collisions(i) = Collisions(i) + dt;          %sum the time between collions per electron
       
       Vt = sqrt(Vx(i)^2 + Vy(i)^2);                %As we loop to check bounds
       tempSum = tempSum + (Mass*Vt^2)/(2*C.kb);    %we might aswell do the temp
                                                    %cacluations
       
       if i <= numVisable
           figure(1);
           subplot(2,1,1);
           plot([xp(i) x(i)], [yp(i) y(i)],'color',colorVec(i,:));
       end
    end
   
    avgTemp = tempSum/numElect;         %evaluate the avg temp of the system
    Temp = [prevTemp avgTemp];          %takes two points to make a line
    Time = [(t-dt) t];                  %the previous temp and the previous
    figure(1);                          %time should line up, so t-dt is the
    subplot(2,1,2);                     %previous temp                                       
    plot(Time, Temp, '-', 'color', colorVec(1,:));
    
    %fprintf('time: %g (%5.2g %%)\n', t, t/TStop*100);
    
    prevTemp = avgTemp;                 
    avgTemp = 0;
    tempSum = 0;
    pause(0.001);
    t = t + dt;
end

%mean free path calculation
%for i = 1:numElect
%    xi(i) = x(i) - xi(i);
%    yi(i) = y(i) - yi(i);   
%end
%avgx = sum(xi - x)/numElect;
%avgy = sum(yi - y)/numElect;
avgx = sum(Vxi - Vx);
avgy = sum(Vyi - Vy);
AvgDist = sum(sqrt(avgDistX.^2 + avgDistY.^2))/numElect;
avgTot = sqrt(avgx^2 + avgy^2)/sqrt(2)*pi*numElect*(AvgDist)^2;

%mean time between collisions
avgCollTime = sumCollTime / numColl;

fprintf('Mean Free Path Calcuated: %g Avg. time between collisions: %g\n', avgTot, avgCollTime); 

if SavePics
    figure(1);
    saveas(gcf, 'ElectronsInSiliconQ2.jpg');

    figure(2);
    saveas(gcf, 'VelocityHistQ2.jpg');
end