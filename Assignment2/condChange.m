C1 = logspace(-2, 1,20);
Currents = zeros(length(C1),1);
extras =1;
for i=1:length(C1)
   voltageFeild2(100, 100, 1, 50, 50, C1(i), 1);
   if extras 
       figure(24);
       subplot(4,5,i);
       mesh(x,y,V2);

       figure(25);
       subplot(4,5,i);
       mesh(x,y,C);
   
       figure(26);
       subplot(4,5,i);
       mesh(x,y,Ex);
       figure(27);
       subplot(4,5,i);
       mesh(x,y,Ey);

       figure(28);
       subplot(4,5,i);
       quiver(x,y,Ex',Ey');
        
       figure(29);
       subplot(4,5,i);
       quiver(x,y,Jx',Jy');
   end

   Currents(i) = curr;
end
condC = Currents;
figure(30);
semilogx(C1,condC);

avgCC = sum(condC)/length(C1);
fprintf('the average current due to changing conductance: %g\n',avgCC);