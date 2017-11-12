clc;
clear all;
close all;
syms x1 x2

%definition of the function
func = 100*(x2-(x1)^2)^2 + (1-x1)^2;
f = matlabFunction(func);

%gradient of the function at xk
gradF = gradient(func);
gradF_eval = matlabFunction(gradF);

%hessian of the function at xk
hessF = hessian(func);
hessF_eval = matlabFunction(hessF);

%setting initial values of backtracking algorithm
threshold = 100;
rho = 0.5 ; %contraction factor
c = 10^(-4); 
x = zeros(2,threshold);
x(:,1) = [-1.2;1];
alpha = ones(1,threshold);
target_error = zeros(1,threshold);
%backtracking algorithm
iterator = 0;

while(iterator <= 100)
    iterator  = iterator + 1
    %evaluating gradient and hessian values at xk
    gradF_value = gradF_eval(x(1,iterator),x(2,iterator));
    hessF_value = hessF_eval(x(1,iterator),x(2,iterator));
    
    %evaluating the search direction pk
    pk = - hessF_value \ gradF_value;
    %evaluating the next point x(k+1)
    x(:,iterator+1) = x(:,iterator)+ alpha(iterator)*pk;
    
    %calculating Armijo's condition
    left_condition = f(x(1,iterator+1),x(2,iterator+1));
    right_condition = f(x(1,iterator),x(2,iterator)) + c*alpha(iterator)*pk'*gradF_value;
 
    %backtracking algorithm
    while(left_condition > right_condition)
         alpha(iterator) = rho*alpha(iterator);
         x(:,iterator+1) = x(:,iterator) + alpha(iterator)*pk;
         left_condition = f(x(1,iterator+1),x(2,iterator+1));
         right_condition = f(x(1,iterator),x(2,iterator)) + c*alpha(iterator)*pk'*gradF_value;
    end
    %condition to check the error rate
    %plotting the graph
   
    target_error(iterator) = abs(abs(f(x(1,iterator+1),x(2,iterator+1))) - abs(f(x(1,iterator),x(2,iterator))));
    if  target_error(iterator) < 10^(-6)
        break;
    end
    title('Error Rate at each iterations');
    semilogy(target_error(1:iterator),'b-o')
    hold on
    grid on;
    str = '| F(x_{k+1})- F(x_k) |';
    text('Interpreter','latex')
    xlabel('Number of Iterations');
    ylabel(str);
end