clearvars

L = 60;                     % Y-direction
W = (2*L/3);                % X-direction

G = sparse(L*W,L*W);

V0 = 1;                     % 1V
diff = 1;

R = ones(W,L);
B = zeros(1,W*L);

%G matrix
for j = 1:L     
    for i = 1:W
        n = i + (j-1)*W;    % I want this to march top left to right, top to bottom
        
        if i == 1           % left side
            G(n,:) = 0;
            G(n,n) = V0;
            B(n) = V0;
            
        elseif i == W       % right side
            G(n,:) = 0;
            B(n) = 0;
            %G(n,n) = V0;
            
        elseif j == 1       % top
            nxp = i + (j-2)*W;
            nxn = i + j*W;
            nyn = i + 1 + (j-1)*W;
            
            G(n,:) = 0;
            %G(n,n) = V0;
            B(n) = (G(n,nyn) - G(n,n))/2.0;
        elseif j == L       % bottom
            nxp = i + (j-2)*W;
            nxn = i + j*W;
            nyp = i - 1 + (j-1)*W;
            
            G(n,:) = 0;
            %G(n,n) = V0;
            B(n) = (G(n,nyp) - G(n,n))/2.0;
        else                % interrior
            nxp = i + (j-2)*W;
            nxn = i + j*W;
            nyp = i - 1 + (j-1)*W;
            nyn = i + 1 + (j-1)*W;
            
            rxn = (R(i,j) + R(i+1))/2.0;
            rxp = (R(i,j) + R(i-1,j))/2.0;
            ryn = (R(i,j) + R(i,j+1))/2.0;
            ryp = (R(i,j) + R(i,j-1))/2.0;
            
            G(n,n) = -(rxn+rxp+ryn+ryp);
            G(n,nxp) = rxp;
            G(n,nxn) = rxn;
            G(n,nyn) = ryn;
            G(n,nyp) = ryp;
        end
    end
end
figure(1);
spy(G);

V = G\B';

