function [ parts ] = electronModling( numElectrons )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Specifing a minimum of one electon to be evaluated
if numElectrons == 0
    return
end

T = 300;
Mo = 9.109 * 10 ^-31;
Mn = 0.26*Mo;
step = 2;
k = 1.381 * 10 ^-23;
vth = sqrt(3*k*T/Mn);

%each Particle is a representation of each electrons starting position
Particle = x:(0-200)*rand() 
           y:(0-100)*rand()
           xVel:(vth)*rand()
           yVel:(vth)*rand()
           
parts = {1,numElectrons};

for p = 1:numElectrons
    %update particle positions, and velocities
     

end


end

