clc;
clear all;
close all;
syms x1 x2

%definition of the function
func = 100*(x2-(x1)^2)^2 + (1-x1)^2;
f = matlabFunction(func);
xk = [-1.2 1.2];

%gradient of the function at xk
gradF = gradient(func);
gradF_eval = matlabFunction(gradF);

%hessian of the function at xk
hessF = hessian(func);
hessF_eval = matlabFunction(hessF);

%setting initial values of backtracking algorithm
threshold = 10000;
rho = 0.5 ; %contraction factor
c = 10^(-4); 
x = zeros(2,threshold);
x(:,1) = [-1.2;1];
alpha = ones(1,threshold);
target_error = zeros(1,threshold);
%backtracking algorithm
iterator = 0;
    title('Descending points each iterations');
    plot(x(1,iterator+1),x(2,iterator+1),'rx')
    hold on
    grid on;
    str = 'X2 values';
    text('Interpreter','latex')
    xlabel('X1 Values');
    ylabel(str);
    xlim([-2,2]);
    ylim([-2,2]);
while(iterator <= 10000)
    iterator  = iterator + 1
    %evaluating gradient and hessian values at xk
    gradF_value = gradF_eval(x(1,iterator),x(2,iterator));
    hessF_value = hessF_eval(x(1,iterator),x(2,iterator));
    
    %evaluating the search direction pk
    pk = - gradF_value /  norm(gradF_value);
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
    
end

title('Descending points each iterations');
    plot(x(1,1:iterator),x(2,1:iterator),'b->')
    hold on
    grid on;
    str = 'X2 values';
    text('Interpreter','latex')
    xlabel('X1 Values');
    ylabel(str);
    xlim([-2,2]);
    ylim([-2,2]);
    legend('Start point', 'Descending points')
    hold off