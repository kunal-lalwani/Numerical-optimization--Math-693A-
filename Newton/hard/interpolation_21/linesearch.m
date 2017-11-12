function [X,min] = linesearch(X)
clc;
close all;
syms x1 x2 

func = (100.*((x2-(x1.^2)).^2))+((1-x1).^2);

f = @(x) (100.*((x(2)-(x(1).^2)).^2))+((1-x(1)).^2);
gradF = @(x) [(100.*((4.*(x(1).^3))-(4.*x(1).*x(2))))+((2.*x(1))-2) ,...
            (100.*((2.*x(2))-(2.*(x(1).^2))))];
        
%Evaluate Phi and Phi' value
phi = @(X,P,alpha) f(X+(alpha.*P)); 
phiDash = @(X,P,alpha) P*gradF(X+(alpha.*P))';  

%function for finding Search Direction
searchDirectionSD = @(x) -x/norm(x);
hessF = @(x) [(100*((12*(x(1)^2))-(4*x(2)))+2) , (100*(-4*x(1))) ; (100*(-4*x(1))) , 200];
searchDirectionNewton = @(x) (- gradF(x) * pinv(hessF(x)));

alpha = [];
storeAlpha = zeros(1,10000);
X = [-1.2 1];      % intitial value
iterator = 0;              % Counter for number of iterations
% pk = searchDirectionSD(gradF(X));      % claculate the steepest direction
pk = searchDirectionNewton(X);             % Calculate the newton direction

points = (X);
alpha_star = 0;
while (norm(gradF(X)) > 10^-6)   % Loop until the gradient is 
                                    % nearly equal to Zero
    initialAlpha = 1;
    X1 = X + (initialAlpha.* pk);
    alpha(1) = 0;
   % Initial alpha for Steepest Direction
%     alpha(2) = (pk * gradF(X)')/(searchDirectionSD(gradF(X1)) * gradF(X1)');
    
   % Initial alpha for Newton
   alpha(2)=1; 

    alphaMax = 1;
    c1 = 10^-4; 
    c2 = 0.9;
    i =2;
    while (true)
        f_eval = phi(X,pk,alpha(i));
        % Armijo's condition check
        if or((f_eval>(phi(X,pk,0)+(c1.*alpha(i).*phiDash(X,pk,0)))),and(i>2,(f_eval>=phi(X,pk,alpha(i-1))))) 
            %Taking a larger step
            alpha_star = czoom(X,pk,alpha(i-1),alpha(i)); 
            break;
        end
        grad_eval = phiDash(X,pk,alpha(i));
        
        %Curvature Wolfe condition check
        if norm(grad_eval)<= (-c2*phiDash(X,pk,0)) 
            alpha_star = alpha(i);
            break;
        end
        
        %Check if Going uphill direction
        if phiDash(X,pk,alpha(i))>=0     
            alpha_star = czoom(X,pk,alpha(i),alpha(i-1));
            break;
        end
        alpha(i+1)= (alphaMax-alpha(i))/2;
        i = i+1;
    end
    storeAlpha(iterator+1) = alpha_star;
    % updating the initial value to point to new direction
    X = X + (alpha_star .* pk);      
    
    % Calculate the New SearchDirection for Steepest Descent
%     pk = searchDirectionSD(gradF(X));                              
    pk = searchDirectionNewton(X);     % Calculate the NEW newton direction
    
    iterator = iterator + 1;
    %adding the point traversed in the Matrix
    points = [points;X];
    
    %displaying the information for reference
    fprintf('Iteration number : %d\n',iterator);
%     fprintf('-----------------------\n');
    fprintf('Value of the variables : [%f %f]\n',X');
%     fprintf('Value of the function : %f\n',f(X));
%     fprintf('Value of the gradient : [%f %f]\n', gradF(X)');
%     fprintf('Value of Alpha : %f\n\n\n',alpha_star);
       
end
    %Plotting Step length Values at each iterations
    figure
    loglog(storeAlpha,'b-o')
    hold on
    grid on;
    str = 'Step length \alpha';
    text('Interpreter','latex')
    xlabel('Number of Iterations');
    ylabel(str);
    title('Step length Values at each iterations');
    
    % Plotting contours and descent path.
    figure
    points = [points;X];          
    [x, y] = meshgrid(-1.5:0.05:1.5 , -1.3:0.05:1.5);
    z = (100.*((y-(x.^2)).^2))+((1-x).^2);
    contour(x,y,z,50);
    hold on;
    plot(points(:,1),points(:,2),'r-','LineWidth',1.5);
    plot(points(end),'bx','LineWidth',1.5,'MarkerSize',10);
    hold off;
    legend('Rosenbrock Function','Descending Point','Optimum Point')
    title(['f(x,y) = ' char(func)]);
    xlabel('x_{1}');
    ylabel('x_{2}');
    min = f(X); 
    grid on
    
end       
    