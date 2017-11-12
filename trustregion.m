% MODIFIED ROSENBROCK FUNCTION
% ############################
% 
% f(x1,x2) = (10*((x(2)-(x(1)^2))^2))+((1-x(1))^2)
% 
%
% Functions Used : 
% ###############
% 
% 1. hessf(X)   : Hessian at X. 
%                      
% 2. df   : Gradient at X.
%                     
% 3. f(X)       : Fuction evaluation at the given point X.
%
% 4. pFullStep(X)    : Full step from the point X. 
%                        
% 5. pSteepestDirection(X)     : Step towards Steepest Descent Direction from X.



function [pDogleg,pOptimal] = trustregion()
clc;
close all;
% X = [0 , 0.5];
X = [0 , -1];

%definition of function, it's gradient and hessian
f = @(x)(10.*((x(2)-(x(1).^2)).^2))+((1-x(1)).^2);
df = @(x) [ 2*x(1) - 40*x(1)*(- x(1)^2 + x(2)) - 2;...
                - 20*x(1)^2 + 20*x(2)];
hessf = @(x) [120*x(1)^2 - 40*x(2) + 2, -40*x(1) ; -40*x(1), 20];

pFullStep = @(x) -hessf(x)\df(x);
pSteepestDirection = @(x) -df(x)/ norm(df(x));

%------------Variable Initialisation---------------------------------------
tau = 0;        
pDogleg = [];      
pOptimal = [];      % Array for Optimal path
delta = 0.1:2/10:1.98;      % 10 values between 0 to 2
L = [];
q = [];
pDoglegArray = [];        %Array for storing dogleg values

%----------------- Dogleg Algorithm----------------------------------------
for i=1:length(delta)           
    if(norm(pSteepestDirection(X))>=delta(i))
    %Steepest Descent Direction
        pDogleg = delta(i) .* (pSteepestDirection(X)/norm(pSteepestDirection(X)));   
        pDoglegArray = [pDoglegArray;pDogleg'];
   
    elseif(norm(pFullStep(X))<=delta(i))
    %Full step (newton) Direction
        pDogleg = pFullStep(X);     
        pDoglegArray = [pDoglegArray;pDogleg'];
    
    else
        %Dogleg Path
        t = vpasolve(norm(pSteepestDirection(X)+(tau-1).*(pFullStep(X)-pSteepestDirection(X))).^2 == delta.^2,tau,[1,2]);
        pDogleg(i) = pSteepestDirection(X)+(t-1).*(pFullStep(X)-pSteepestDirection(X));   
        pDoglegArray = [pDoglegArray;pDogleg'];
        
    end  
    
end

%------- grid for plotting model contour-----------------------------------
x = -0.8:0.05:1.1;      
y = -0.8:0.05:1.1;
z = [];

for i=1:length(x)
    
    for j=1:length(y)
        %Quadratic model
        b = [x(i),y(j)];
        z(i,j) = f(X)+b*df(X)+((b*hessf(X)*b').*0.5);  
        
    end
    
end


for j=1:length(delta)   
    lambda = [];
    lambda(1) = 0;
    k = 1;
    p = [];
    %-----------Exact trust region Algorithm------------------------------
    while(1)        
        while(1)
            %Cholesky Decomposition
            [L, flag] = chol(hessf(X)+(lambda(k).*eye(2))); 
            if(flag)
                lambda(k) = lambda(k) + 0.02;
            else
                break;
            end
        end
        
     p = [p ; (-(L'*L)\df(X))']; 
     q = [q ; (L'\p(k,:)')'];
     
     lambda(k+1) = lambda(k) + ((norm(p(k,:))/norm(q(k,:))).^2)...
                                .*((norm(p(k,:))-delta(j))/delta(j)); 
     %termination condition
     if((lambda(k+1) - lambda(k)) < 10^-4)  
         break;
     end
     
     k = k+1;
    
    end
    
    pOptimal = [pOptimal ;p(k,:)];

end

%-------------------plot contours------------------------------------------
figure(001)     
contour(x,y,z',30);
hold on;
s1 = pSteepestDirection(X);
s2 = pFullStep(X);
plot([0,s1(1)],[0,s1(2)],'-rx');
plot([0,s2(1)],[0,s2(2)],'-bo');

pDoglegArray = [[0,0];pDoglegArray];
pOptimal = [[0,0];pOptimal];
pOptimal = [pOptimal ; pFullStep(X)'];

plot(pDoglegArray(:,1),pDoglegArray(:,2),'-*k','LineWidth',1.1);

plot(pOptimal(:,1),pOptimal(:,2),'--m','LineWidth',1.1);
grid on;
hold off;

title('Dogleg path vs Optimal Path');
legend('Model','Steepest Descent Direction','Newton Direction','Dogleg Path','Optimal Path',...
       'Location','NorthEastOutside')
pDogleg = pDoglegArray;
colorbar;
% xlim([-0.4 0.4]);
% ylim([-1.5 0.2]);
xlim([-0.1 0.15]);
ylim([-0.2 1.2]);
end