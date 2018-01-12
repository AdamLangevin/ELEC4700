function [ vth  ] = electronModling( numElectrons,  Mo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Specifing a minimum of one electon to be evaluated
if numElectrons == 0
    vth = 0;
    return
end

T = 300;
Mn = 0.26*Mo;
%each Particle is a representation of each electrons starting position
Particle = x:(0-200)*rand() 
           y:(0-100)*rand()
           xVel:(



end

