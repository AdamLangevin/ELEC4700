function [ parts ] = electronModling( numElectrons, T )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Specifing a minimum of one electon to be evaluated
if numElectrons == 0
    return
end

Mo = 9.109 * 10 ^-31;
Mn = 0.26*Mo;
k = 1.381 * 10 ^-23;
vth = sqrt(3*k*T/Mn);

taumn = 0.2*e-12;

L = 200*e-9;
W = 100*e-9;

xp(1,:) = linspace(0, L, numElectrons);
yp(1,:) = linspace(0, W, numElectrons);

x(numElectrons + 1:numElectrons) = xp - L/2;
y(numElectrons + 1:numElectrons) = y(1) - W/2;

% Each Particle is a representation of each electrons starting position
% this is the important information, but this structure is difficult to use
%Particle = x:(0-200)*rand() 
%           y:(0-100)*rand()
%           xVel:(vth)*rand()
%           yVel:(vth)*rand()
           
parts = {1,numElectrons};

for i = 1:numElectrons
    %initilize particle positions, and velocities
    

end


end

