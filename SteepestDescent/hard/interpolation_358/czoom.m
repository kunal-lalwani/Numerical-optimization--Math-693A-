function alpha = czoom(X,P,alpha_low,alpha_high)

f = @(x)(100.*((x(2)-(x(1).^2)).^2))+((1-x(1)).^2);
f1 = @(x) [(100.*((4.*(x(1).^3))-(4.*x(1).*x(2))))+((2.*x(1))-2) , ...
           (100.*((2.*x(2))-(2.*(x(1).^2))))];

% calculate Phi and Phi'
phi = @(X,P,alpha) f(X+(alpha.*P)); 
phiDash = @(X,P,alpha) P*f1(X+(alpha.*P))';

c1 = 10^-4;
c2 = 0.9;
while true     
    %functions for calculating phi and phi' at different alphas
    phi_low = phi(X,P,alpha_low);
    phi_high = phi(X,P,alpha_high);
    phiDash_low = phiDash(X,P,alpha_low);
    phiDash_high = phiDash(X,P,alpha_high);
    phi0 = phi(X,P,0);
    phiDash0 = phiDash(X,P,0);
    
    %Calculate alpha by interpolation
    alphaj = interpolations(alpha_low,alpha_high,phi_low,phi_high,phiDash_low,phiDash_high);
%     alphaj = (alpha_low+alpha_high)/2;  
    %Check if alpha satisfies Armijo's condition
    a1 = phi(X,P,alphaj);
    if or((a1 > (phi0+(c1.*alphaj.*phiDash0))),(a1 >= phi_low))
        alpha_high = alphaj;
    else
        a_1 = phiDash(X,P,alphaj);
        %Strong Wolfe's curvature condition
        if norm(a_1)<= (-c2*phiDash0)
            alpha = alphaj;
            return;
        end
        if a_1*(alpha_high-alpha_low) >= 0
            alpha_high = alpha_low;
        end
        alpha_low = alphaj;
    end
end

        
