function  PA3(N, M, r, l, t, b, insulTB, insulRL)

n = N;
m = M;

right = r;
left = l;
top = t;
bot = b;

if insulTB == 1
    
end

if insulRL ==1 
    
end

global V;

for i = 1:m
    for j = 1:n
        if i-1 <= 0                             %top wall
            if insulTB == 0
                if j-1 <= 0
                    V(i,j) = (V(i+1,j) + top + V(i,j+1) + left)/4;

                elseif j+1 > n
                    V(i,j) = (V(i+1,j) + right + top + V(i,j-1))/4;

                else
                    V(i,j) = (V(i+1,j) + top + V(i,j+1) + V(i,j-1))/4;

                end
            else
               if j-1 <= 0
                    V(i,j) = (V(i+1,j) + V(i,j) + V(i,j+1) + left)/4;

                elseif j+1 > n
                    V(i,j) = (V(i+1,j) + right + V(i,j) + V(i,j-1))/4;

                else
                    V(i,j) = (V(i+1,j) + V(i,j) + V(i,j+1) + V(i,j-1))/4;

                end 
            end
        end
        if i+1 > m                              %bottom wall
            if insulTB == 0
                if j-1 <= 0
                    V(i,j) = (V(i-1,j) + bot + V(i,j+1) + left)/4;

                elseif j+1 > n 
                    V(i,j) = (bot + V(i-1,j) + right + V(i,j-1))/4;

                else
                    V(i,j) = (bot + V(i-1,j) + V(i,j+1) + V(i,j-1))/4;

                end
            else
                if j-1 <= 0
                    V(i,j) = (V(i-1,j) + V(i,j) + V(i,j+1) + left)/4;

                elseif j+1 > n 
                    V(i,j) = (V(i,j) + V(i-1,j) + right + V(i,j-1))/4;

                else
                    V(i,j) = (V(i,j) + V(i-1,j) + V(i,j+1) + V(i,j-1))/4;

                end
            end
        end
        if j-1 <= 0                            %left nodes
            if insulTB ==0
                if i-1<=0
                    %%V(i,j) = (V(i,j) + V(i+1,j) + V(i,j+1) + left)/4;

                elseif i+1>m
                    %V(i,j) = (V(i-1,j) + V(i,j) + V(i,j+1) + left)/4;

                else
                    V(i,j) = (V(i-1,j) + V(i+1,j) + V(i,j+1) + left)/4;

                end 
            else
                if i-1<=0
                    %V(i,j) = (top + V(i+1,j) + V(i,j+1) + left)/4;

                elseif i+1>m
                    %V(i,j) = (V(i-1,j) + bot + V(i,j+1) + left)/4;

                else
                    V(i,j) = (V(i-1,j) + V(i+1,j) + V(i,j+1) + left)/4;

                end  
            end
        end
        if j+1>n                               %right nodes
            if insulTB ==0
                 if i-1<=0
                    %V(i,j) = (V(i,j) + V(i+1,j) + right + V(i,j-1))/4;

                elseif i+1>m
                    %V(i,j) = (V(i-1,j) + V(i,j) + right + V(i,j-1))/4;

                else
                    V(i,j) = (V(i-1,j) + V(i+1,j) + right + V(i,j-1))/4;

                end
            else
                if i-1<=0
                    %V(i,j) = (top + V(i+1,j) + right + V(i,j-1))/4;

                elseif i+1>m
                    %V(i,j) = (V(i-1,j) + bot + right + V(i,j-1))/4;

                else
                    V(i,j) = (V(i-1,j) + V(i+1,j) + right + V(i,j-1))/4;

                end
            end
            
        end                      
        if i>1 && i<n && j>1 && j<n           %interrior nodes
            V(i,j) = (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1))/4;
        end
                
    end
end
end