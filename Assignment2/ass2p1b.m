clearvars

nx = 100;
ny = nx*3/2;
a = 1;
b = 1.5;
numSims = 21;
V0 = 1;
V = zeros(nx,ny);
Vp = zeros(nx,ny);
V(nx,:) = V0;
V(1,:) = V0;
X = linspace(-b/2,b/2,nx);
Y = linspace(0,a,ny);

for n = 1:2:numSims 
    for i = 1:nx
        for j = 1:ny
            if X(i) == -b/2
                V(i,j) = 1;
            elseif X(i) == b/2
                V(i,j) = 1;
            elseif Y(j) == 0
                V(i,j) = 0;
            elseif Y(j) == a
                V(i,j) = 0;
            else
                V(i,j) = (1/n) * (cosh(n*pi*X(i)/a) / cosh(n*pi*b/a)) * sin(n*pi*Y(j)/a) + Vp(i,j);
            end
        end 
    end
    V = (4*V0/pi).*V;
    Vp = V;
    pause(0.0000001);
    figure(1);
    mesh(Y,X,V);
end
