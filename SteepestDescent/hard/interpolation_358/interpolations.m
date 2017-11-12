function alpha = interpolations(alpha_low,alpha_high,phi_low,phi_high,phidash_low,phidash_high)

% Interpolate using Hermite polinomial and find the step length alpha.
d1 = phidash_low+phidash_high-3.*((phi_low-phi_high)/(alpha_low-alpha_high));
d2 = real(sign(alpha_high-alpha_low)*((d1.^2 - (phidash_low.*phidash_high))  ^0.5));
alpha=alpha_high-((alpha_high-alpha_low).*((phidash_high+d2-d1)/(phidash_high-phidash_low+(2.*d2))));

return