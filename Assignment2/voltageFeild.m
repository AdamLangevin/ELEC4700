function [] = voltageFeild (length, width, initialBC, twoEnded)

global V2;
global x;
global y;

Lx = 1;
Ly = 1.5;
nx = width;
ny = length;
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
diff = (1/width)^2;
V0 = initialBC;
V2 = zeros(nx,ny);

G = sparse(nx*ny,nx*ny);
B = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = i + (j-1)*nx;
        
        if i == 1                   %left
           G(n,n) = 1;
           B(1,n) = V0;
           
        elseif i ==nx               %right
           G(n,n) = 1;
           if twoEnded == 1
               B(1,n) = V0;
           else
               B(1,n) = 0;
           end
           
        elseif j == 1               %bottom
            nxn = i + 1 + (j-1)*nx;
            nxp = i - 1 + (j-1)*nx;
            nyn = i + (j)*nx;
            
            G(n,n) = 1;
            
        elseif j == ny              %top side
            nxn = i + 1 + (j-1)*nx;
            nxp = i - 1 + (j-1)*nx;
            nyp = i + (j-2)*nx;
            
            G(n,n) = 1;
            
        else                        %interrior            
            nxn = i + 1 + (j-1)*nx;
            nxp = i - 1 + (j-1)*nx;
            nyp = i + (j-2)*nx;
            nyn = i + (j)*nx;
            
            G(n,n) = -4;
            G(n,nxn) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
            G(n,nyn) = 1;         
        end
    end    
end
%figure(2);
%spy(G);

V = G\B';

for i = 1:nx
    for j = 1:ny
        n = i + (j-1)*nx;
        V2(i,j) = V(n);        
    end
end

%x = linspace(1,nx,nx);
%y = linspace(1,ny,ny);

%mesh(y, x, V2)
end
