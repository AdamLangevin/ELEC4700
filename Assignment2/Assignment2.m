%% Assignment 2
%% Question 1
% Adam Langevin 100935879.
%
% a) At 80 points for the total length, with the boundary conditions V=V0 
% at x=0, and V = 0 at x=L, with the sides being isolated.
%
% b) The case for V=V0 at x=0 and x=L wuth a density of 100 points per unit 
% lenght, the V(x) plot of the voltage as a function of position is in 
% figure 3 and 4.
%
% The size of the meshing of the area shows more detail of the area as the
% number of points calcualted increasses, but the time needed to calculate
% it increases at a much faster rate when in two dimensions than when in
% one. Depending on the problem the mesh size only shows so much, and
% adding more points doesnt add more information. 
%
% The analytical solution provides some insight on the interior of the
% feild, and shows the general shape. Around the sides where the boundary 
% conditions are the functional solution we have created shows better 
% information on how the structure behaves.

%% Question 1 code

global x y V2;

% a)
figure(1)
voltageFeild(80,80,1,0);
plot(x,V2);
title('1-D potential at 80 points');
xlabel('x-position');
ylabel('potential');

% The complimentary surface plot for the same mesh size is,
figure(2)
mesh(y, x, V2);
title('2-D potential at 80 points');
xlabel('y');
ylabel('x');

% b)
figure(3);
voltageFeild(100,100,1,1);
plot(x,V2);
title('The 1-D potential plot in the region');
xlabel('x');
ylabel('Voltage');

% The meshed voltage in the x and y directions.
figure(4);
mesh(y,x,V2);
title(['2-D mesh plots of potential at ', num2str(100), ' points']);
xlabel('y');
ylabel('x');
zlabel('Potential as a function of x and y');

%% Question 2
% a) The Voltage feild V(x,y), the conductivity map, the electric feild, and
% the curent densisty map are in figures 5-10 respecivly.
%
% b) calculating the density of the meshingvs the current density passing
% through each end of the feild. By increasing the amount of points per
% unit the details showed more, creating sharper peaks. The plot of mesh 
% density vs current the current generally increases.
%
% c) Changing only the width of the bottle neck did not significantly 
% change the current through the feild. By turing on the extras, figure 22
% would show that the current density around the ends of the bottle neck
% would be different, but the current it's self would not change.
%
% d) To investigate the changing of the current while the value of the
% conductivity inside the bottle neck "boxes" changes I used a log scale to
% change the conductivity from very low values of 0.01 to 10. The
% resulting graph, figure 30, shows a steady average value around 6.75 fA 
% of current was formed through the material, with rough general increasing
% trend.

%% Question 2 Code

% to plot each graph set extras to 1.
% a)
extras = 0;
global Ex Ey Jx Jy C;
voltageFeild2(100,100,1,40,30,0.01,1);

figure(5);
mesh(x,y,V2);
title('Potential in a plane');
xlabel('x');
ylabel('Y');
zlabel('Potential');

figure(6);
mesh(x,y,C);
title('Conductivity Map');
xlabel('x');
ylabel('Y');
zlabel('Conductivity');

figure(7);
mesh(x,y,Ex);
title('Electric feild in the X direction');
xlabel('x');
ylabel('Y');
zlabel('Electrical Potential Energy');

figure(8);
mesh(x,y,Ey);
title('Electric feild in the Y direction');
xlabel('x');
ylabel('Y');
zlabel('Electrical Potential Energy');

figure(9);
quiver(x,y,Ex',Ey');
title('Electric feild in the region');
xlabel('x');
ylabel('Y');

figure(10);
quiver(x,y,Jx',Jy');
title('Current Density in the region');
xlabel('x');
ylabel('Y');

% b)
X = [10 20 30 40 50 75 100 200];
Currents = zeros(length(X),1);

for i = 1:length(X)
   voltageFeild2(X(i), X(i), 1, 0.2*X(i), 0.2*X(i), 0.01, 1);
   
   if extras
       figure(11);
       subplot(4,4,i);
       mesh(x,y,V2);
       title('Potential in a plane');
       xlabel('x');
       ylabel('Y');
       zlabel('Potential');
   
       figure(12);
       subplot(3,3,i);
       mesh(x,y,C);
       title('Conductivity Map');
       xlabel('x');
       ylabel('Y');
       zlabel('Conductivity');

       figure(13);
       subplot(3,3,i);
       mesh(x,y,Ex);
       title('Electric feild in the X direction');
       xlabel('x');
       ylabel('Y');
       zlabel('Electrical Potential Energy');
       
       figure(14);
       subplot(3,3,i);
       mesh(x,y,Ey);
       title('Electric feild in the Y direction');
       xlabel('x');
       ylabel('Y');
       zlabel('Electrical Potential Energy');
       
       figure(15);
       subplot(3,3,i);
       quiver(x,y,Ex',Ey');
       title('Electric feild in the region');
       xlabel('x');
       ylabel('Y');
        
       figure(16);
       subplot(3,3,i);
       quiver(x,y,Jx',Jy');
       title('Current Density in the region');
       xlabel('x');
       ylabel('Y');
   end
   
   Currents(i) = curr;
end
meshC = Currents;
figure(17);
plot(X,meshC);
title('Change in Current as Mesh Density Increases');
xlabel('Number of Mesh points');
ylabel('Difference of Current Between Ends');

avgMC = sum(meshC)/length(meshC);
fprintf('the average current due to the changing mesh size: %g\n',avgMC);

% c)
nx = 100;
W = [nx*0.1 nx*0.15 nx*0.2 nx*0.25 nx*0.3 nx*0.35 nx*0.4 nx*0.45 nx*0.5 ...
    nx*0.55 nx*0.6 nx*0.65 nx*0.7 nx*0.75 nx*0.8 nx*0.85 nx*0.9 nx*0.95];
Currents = zeros(length(W),1);

for i = 1:length(W)
   voltageFeild2(nx, nx, 1, W(i), 40, 0.01, 1);
   if extras
       figure(18);
       subplot(4,5,i);
       mesh(x,y,V2);
       title('Potential in a plane');
       xlabel('x');
       ylabel('Y');
       zlabel('Potential');

       figure(19);
       subplot(4,5,i);
       mesh(x,y,C);
       title('Conductivity Map');
       xlabel('x');
       ylabel('Y');
       zlabel('Conductivity');

       figure(20);
       subplot(4,5,i);
       mesh(x,y,Ex);
       title('Electric feild in the X direction');
       xlabel('x');
       ylabel('Y');
       zlabel('Electrical Potential Energy');

       figure(10);
       subplot(4,5,i);
       mesh(x,y,Ey);
       title('Electric feild in the Y direction');
       xlabel('x');
       ylabel('Y');
       zlabel('Electrical Potential Energy');
       
       figure(21);
       subplot(4,5,i);
       quiver(x,y,Ex',Ey');
       title('Electric feild in the region');
       xlabel('x');
       ylabel('Y');
        
       figure(22);
       subplot(4,5,i);
       quiver(x,y,Jx',Jy');
       title('Current Density in the region');
       xlabel('x');
       ylabel('Y');
   end
  
   Currents(i) = curr;
end
widthC = Currents;
figure(23);
plot(W,widthC);
title('Change in Current While Resistive Feature Changes in Size');
xlabel('Change in Width of Apature');
ylabel('Current');

avgWC = sum(widthC)/length(W);
fprintf('the average current for the chaning width: %g\n', avgWC);

% d)
C1 = logspace(-2, 1,20);
Currents = zeros(length(C1),1);

for i=1:length(C1)
   voltageFeild2(nx, nx, 1, 50, 50, C1(i), 1);
   if extras 
       figure(24);
       subplot(4,5,i);
       mesh(x,y,V2);
       title('Potential in a plane');
       xlabel('x');
       ylabel('Y');
       zlabel('Potential');

       figure(25);
       subplot(4,5,i);
       mesh(x,y,C);
       title('Conductivity Map');
       xlabel('x');
       ylabel('Y');
       zlabel('Conductivity');

       figure(26);
       subplot(4,5,i);
       mesh(x,y,Ex);
       title('Electric feild in the X direction');
       xlabel('x');
       ylabel('Y');
       zlabel('Electrical Potential Energy');

       figure(27);
       subplot(4,5,i);
       mesh(x,y,Ey);
       title('Electric feild in the Y direction');
       xlabel('x');
       ylabel('Y');
       zlabel('Electrical Potential Energy');
       
       figure(28);
       subplot(4,5,i);
       quiver(x,y,Ex',Ey');
       title('Electric feild in the region');
       xlabel('x');
       ylabel('Y');
        
       figure(29);
       subplot(4,5,i);
       quiver(x,y,Jx',Jy');
       title('Current Density in the region');
       xlabel('x');
       ylabel('Y');
   end

   Currents(i) = curr;
end
condC = Currents;
figure(30);
semilogx(C1,condC);
title('Current Change in Region While Conductivity Increases');
xlabel('Conductivity in Box Regions');
ylabel('Current');

avgCC = sum(condC)/length(C1);
fprintf('the average current due to changing conductance: %g\n',avgCC);