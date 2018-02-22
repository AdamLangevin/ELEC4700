
clearvars 
clearvars -Global

global C x y V2 Ex Ey Jx Jy curr;

extras = 0; %this is here to suppress some of the plots not really needed to see the chages

% part b
X = [10 20 30 40 50 100 200 300 400 500 1000 1100 1200 1300 1400 1500];
Currents = zeros(length(X),1);

for i = 1:length(X)
   voltageFeild2(X(i), X(i), 1, 0.2*X(i), 0.2*X(i), 0.01, 1);
   figure(1);
   subplot(4,4,i);
   mesh(x,y,V2);
   
   if extras
       figure(2);
       subplot(4,4,i);
       mesh(x,y,C);

       figure(3);
       subplot(4,4,i);
       mesh(x,y,Ex);
       figure(4);
       subplot(4,4,i);
       mesh(x,y,Ey);

       figure(5);
       subplot(4,4,i);
       quiver(x,y,Ex',Ey');

       figure(6);
       subplot(4,4,i);
       quiver(x,y,Jx',Jy');
   end
   
   Currents(i) = curr;
end
meshC = Currents;
figure(19);
plot(X,meshC);

% part c
% changing only width
nx = 100;
W = [nx*0.1 nx*0.15 nx*0.2 nx*0.25 nx*0.3 nx*0.35 nx*0.4 nx*0.45 nx*0.5 nx*0.55 nx*0.6 nx*0.65 nx*0.7 nx*0.75 nx*0.8 nx*0.85 nx*0.9 nx*0.95];
Currents = zeros(length(W),1);

for i = 1:length(W)
   voltageFeild2(nx, nx, 1, W(i), 0.2*nx, 0.01, 1);
   figure(7);
   subplot(4,5,i);
   mesh(x,y,V2);
   
   figure(8);
   subplot(4,5,i);
   mesh(x,y,C);
   
   if extras
       figure(9);
       subplot(4,5,i);
       mesh(x,y,Ex);
       figure(10);
       subplot(4,5,i);
       mesh(x,y,Ey);

       figure(11);
       subplot(4,5,i);
       quiver(x,y,Ex',Ey');
   end
   
   figure(12);
   subplot(4,5,i);
   quiver(x,y,Jx',Jy');
   
   Currents(i) = curr;
end
widthC = Currents;
figure(20);
plot(W,widthC);

% part c
C1 = logspace(-4,0,12);
Currents = zeros(length(C1),1);

for i=1:length(C1)
   voltageFeild2(nx, nx, 1, 20, 20, C1(i), 1);
   figure(13);
   subplot(4,3,i);
   mesh(x,y,V2);
   
   figure(14);
   subplot(4,3,i);
   mesh(x,y,C);
   
   if extras 
       figure(15);
       subplot(4,3,i);
       mesh(x,y,Ex);
       figure(16);
       subplot(4,3,i);
       mesh(x,y,Ey);

       figure(17);
       subplot(4,3,i);
       quiver(x,y,Ex',Ey');
   end

   figure(18);
   subplot(4,3,i);
   quiver(x,y,Jx',Jy');
   
   Currents(i) = curr;
end
condC = Currents;
figure(21);
plot(C1,condC);

%part a
%voltageFeild(nx,nx,1,
