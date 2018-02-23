clearvars

global V2;
global x;
global y;

X = [10 20 30 40 50 60 70 80 90 100 120 150];

% part 1a, for several 2-D topographies.
figure(1);
for i = 1:length(X)
    voltageFeild(X(i),X(i),1,0);
    subplot(4,3,i);
    mesh(y, x, V2);
    title(['2-D potential at ', num2str(X(i)), ' points']);
    xlabel('y');
    ylabel('x');
end  


figure(2);
for i=1:length(X)
    voltageFeild(X(i),X(i),1,0);
    subplot(4,3,i);
    plot(x,V2);
    title(['1-D potential at ',num2str(X(i)), ' points']);
    xlabel('x-position');
    ylabel('potential');
end

figure(3);
voltageFeild(150,150,1,1);
mesh(y,x,V2);
title(['2-D mesh plots of potential at ', num2str(150), ' points']);
xlabel('y');
ylabel('x');
zlabel('Potential as a function of x and y');
