% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = linspace(3.305,3.34,10);          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute

for i=1:length(n1)
[x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1(i),n2,n3],[h1,h2,h3], ...
                                            rh,rw,side,dx,dy); 

% First consider the fundamental TE mode:

[Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A');

fprintf(1,'neff = %.6f\n',neff);

    figure(4);
    subplot(4,5,i);
    s = contourmode(x,y,20*log10(Hx));
    s.EdgeColor = 'none';
    title(['Hx (TE mode ) ', i]); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end

    subplot(4,5,i+10);
    s = contourmode(x,y,20*log10(Hy));
    s.EdgeColor = 'none';
    title(['Hy (TE mode)', i]); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
end

