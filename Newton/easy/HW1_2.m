clc;
syms x1 x2

%definition of the function
f = @(x1,x2) 100*(x2-(x1)^2)^2 + (1-x1)^2;
xk =  [1.2 1.2]'; %starting point

%gradient of the function at xk
gradF = gradient(f, [x1 x2]);
gradF_value = double(subs(gradF, {x1, x2}, xk'));

%hessian of the function at xk
hessF = hessian(f, [x1 x2]);
hessF_value = double(subs(hessF, {x1, x2}, xk'));

%absolute value of the gradient and hessian
% abs_gradF = double(norm(gradF_value));
% abs_hessF = double(norm(hessF_value));

%search direction from newton's method
pk = - inv(hessF_value) * gradF_value; %search direction

%setting initial values of backtracking algorithm
alpha_bar = 1; %step length
rho = 0.5 ; %contraction factor
c = 10^(-4); 
% alpha= alpha_bar;

%backtracking algorithm
iterator = 0;

while(iterator <= 100)
   
    
    alpha = alpha_bar;
    x_prev = xk; %previous value of xk
    x_next = xk + (alpha * pk);%value of xk+1
   
    
    %calculating Armijo's condition
    left_condition = double(subs(f, {x1, x2}, x_next' )) ;
    right_condition = double(subs(f, {x1, x2}, xk')) + c*alpha*pk'*gradF_value;
    
    %backtracking algorithm
    while(left_condition > right_condition)
        alpha = alpha * rho;
        
        %recalculating gradient and hessian
        gradF = gradient(f, [x1 x2]);
        gradF_value = double(subs(gradF, {x1, x2}, xk'));
        hessF = hessian(f, [x1 x2]);
        hessF_value = double(subs(hessF, {x1, x2}, xk'));
        
        pk = - inv(hessF_value) * gradF_value; %recalculating search direction
        right_condition = double(subs(f, {x1, x2}, xk')) + c*alpha*pk'*gradF_value;
        
        x_next = xk + alpha * pk;
        left_condition = double(subs(f, {x1, x2}, x_next'));
        xk = x_next;
    end
    %setting step length of each iteration
    alpha_k = alpha;
    
    func_prev= double(subs(f, [x1, x2], x_prev'));
    func_next= double(subs(f, [x1, x2], x_prev' + alpha * pk'));
    %error limit |xk+1 - xk|
    if_cond = abs((func_next) - (func_prev));
    if( if_cond < 10^(-6))
        break;
    end
%     figure(1.1)
    value = double(subs(f, {x1, x2}, xk'))  
    title('Step Length Values');
    plot(iterator, alpha, 'b-o','MarkerIndices',1:5:length(f))
    ylim([0,1])
    xlim([0 10])
    hold on
    grid;
    str = 'Step length \alpha';
    text('Interpreter','latex')
    xlabel('Number of Iterations');
    ylabel(str);
    
    xk
    xk = x_prev + alpha * pk;
    iterator  = iterator + 1
end