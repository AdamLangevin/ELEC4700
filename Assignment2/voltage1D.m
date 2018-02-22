function [] = voltage1D(width, initialBC)

global V2;
global x;

nx = width;
diff = 1/nx;
V0 = initialBC;

V2 = zeros(nx,1);
G = zeros(nx);

for i = 1:nx
    if i == 1
        G(i,i) = 1;
    elseif i == nx
        G(i,i) = 1;
    else 
        G(i,i) = -2;
        G(i,i+1) = 1;
        G(i+1,1) = 1;
    end  
end

B = zeros(nx,1);
B(1) = 0;
B(nx) = V0;

V2 = G\B;
x = linspace(0,1,nx);

end

