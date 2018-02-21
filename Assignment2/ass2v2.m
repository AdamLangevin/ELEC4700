clearvars

nx = 100;
ny = nx*3/2;
diff = 1;
V0 = 1;

G = sparse(nx*ny,nx*ny);
B = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = i + (j-1)*nx;
        
        if i == 1                   %left
           G(n,:) = 0;
           G(n,n) = 1;
           B(1,n) = V0;
           
        elseif i ==nx               %right
           G(n,:) = 0;
           G(n,n) = 1;
           B(1,n) = V0;
           
        elseif j == 1               %bottom
            nxn = i + (j-2)*nx;
            nxp = i + (j)*nx;
            nyp = i + 1 +(j-1)*nx;
            
            G(:,n) = 0;
            G(n,n) = 1;
            B(1,n) = 0;
            
        elseif j == ny              %top side
            nxn = i + (j-2)*nx;
            nyn = i - 1 + (j-1)*nx;
            nyp = i + 1 +(j-1)*nx; 
            
            G(:,n) = 0;
            G(n,n) = 1;
            B(1,n) = 0;
        else                        %interrior            
            nxn = i + (j-2)*nx;
            nxp = i + (j)*nx;
            nyp = i + 1 +(j-1)*nx;
            nyn = i - 1 + (j-1)*nx;
            
            G(n,n) = -4/diff;
            G(n,nxn) = 1/diff;
            G(n,nxp) = 1/diff;
            G(n,nyp) = 1/diff;
            G(n,nyn) = 1/diff;
            
        end
    end    
end
figure(1);
spy(G);

V = G\B';

for i = 1:nx
    for j = 1:ny
        n = i + (j-1)*nx;
        V2(i,j) = V(n);        
    end
end

x = linspace(1,nx,nx);
y = linspace(1,ny,ny);

mesh(y, x, V2)

