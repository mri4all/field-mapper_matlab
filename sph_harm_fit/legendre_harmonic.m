function out = legendre_harmonic(n, m, pos_x, pos_y, pos_z)
% returns harmonic field for the required solid harmonic addressed by 
% n, m based on the legendre() Matlab function
% Positive m values correspond to cosine component and negative to sine
%
% returned fields will eventually follow RRI's convention
% pos_... can be both value and vector/matrix

r2=pos_x.^2+pos_y.^2+pos_z.^2;
r=r2.^0.5;
phi=atan2(pos_y, pos_x);
cos_theta=cos(atan2((pos_x.^2+pos_y.^2).^0.5, pos_z));
%cos_theta=pos_z./r;

if m>=0,
    c=1;
else
    c=0;
    m=-m;
end

Ymn=legendre(n,cos_theta);

rri_norm=factorial(n+m+1)/factorial(n-m)/ffactorial(2*m);

out=(n+m+1)*r.^(n).*(cos(m*phi)*c+sin(m*phi)*(1-c)).*reshape(Ymn(m+1,:),size(pos_x))/rri_norm;
