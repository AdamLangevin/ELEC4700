function [] = voltageFeild2(length, width, initialBC, apaW, apaL, C1, C2)

global V2;
global x;
global y;
global Ex Ey Jx Jy C B curr;

Lx = 1;
Ly = 1.5;
nx = width;
ny = length;
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);

%diff = (1/width)^2;

V0 = initialBC;
V2 = zeros(nx,ny);

G = sparse(nx*ny,nx*ny);
B = zeros(1,nx*ny);

apatureWidth = apaW;
apatureLength = apaL;

%set up the conductivity map
C = ones(nx,ny);
for i=1:nx
    for j=1:ny
        if j < (nx/2 - apatureWidth/2) && i > (ny/2 - apatureLength/2) && i < (ny/2 + apatureLength/2 + 1)
            C(i,j) = C1; %inside the top box
        elseif j > (nx/2 + apatureWidth/2 + 1) && i<nx+1 && i >(ny/2 - apatureLength/2) && i < (ny/2 + apatureLength/2 + 1)
           C(i,j) = C1; %inside the bottom box
        else
           C(i,j) = C2;
        end
    end
end

%set up the BCs
for i=1:nx
    for j=1:ny
        n = i +(j-1)*nx;
        if i == 1 || i == nx
            B(n) = V0;
        elseif j == 1 || j == ny
            B(n) = 0;
        end           
    end
end

%G matrix
for i = 1:nx
    for j = 1:ny
       n = j + (i-1)*ny;
        
        if i == 1                   %left
           G(n,n) = 1;
           
        elseif i ==nx               %right
           G(n,n) = 1;
             
        elseif j == 1               %bottom
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (C(i,j) + C(i-1,j))*0.5;
            rxp = (C(i,j) + C(i+1,j))*0.5;
            ryp = (C(i,j) + C(i,j+1))*0.5;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;
            
        elseif j == ny              %top side
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            
            rxm = (C(i,j) + C(i-1,j))*0.5;
            rxp = (C(i,j) + C(i+1,j))*0.5;
            rym = (C(i,j) + C(i,j-1))*0.5;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
            
        else                        %interrior            
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (C(i,j) + C(i-1,j))*0.5;
            rxp = (C(i,j) + C(i+1,j))*0.5;
            rym = (C(i,j) + C(i,j-1))*0.5;
            ryp = (C(i,j) + C(i,j+1))*0.5;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;         
        end
    end    
end

V = G\B';

%reassign to an x,y coordinate system
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*nx;
        V2(i,j) = V(n);        
    end
end

Ex = zeros(nx,ny);
Ey = zeros(nx,ny);

%calcualte the electric feilds
for i = 1:nx
 for j = 1:ny
     if i == 1
         Ex(i,j) = (V2(i+1,j) - V2(i,j));
     elseif i == nx
         Ex(i,j) = (V2(i,j) - V2(i-1,j));
     else
         Ex(i,j) = (V2(i+1,j) - V2(i-1, j))/2.0;
     end
     if j == 1
         Ey(i,j) = (V2(i,j+1) - V2(i,j));
     elseif j == ny
         Ey(i,j) = (V2(i,j) - V2(i, j-1));
     else
         Ey(i,j) = (V2(i,j+1) - V2(i,j-1))/2.0;
     end
 end
end

Ex = -Ex;
Ey = -Ey;

%Current density
Jx = C.*Ex;
Jy = C.*Ey;

SC0 = sum(Jx(1,:));
SC1 = sum(Jx(nx,:));
curr = (SC0 + SC1)*0.5;
end

