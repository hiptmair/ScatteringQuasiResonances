% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
function [x,w]=gaussquad(a,b,n)
% n-pt Gauss quadrature nodes and weights for the interval [a,b]
% Relies on Golub-Welsh algorithm
od = zeros(n-1,1);
for i=1:(n-1), od(i)=i/sqrt(4*i*i-1); end
J=diag(od,-1)+diag(od,1); [ev,ew]=eig(J);
for i=1:n, ev(:,i)./norm(ev(:,i)); end
x=diag(ew)'; w=(2*(ev(1,:).*ev(1,:)));
% Affine transformation of points 
x = a+(b-a)*(x+1)/2;
% Transformation of weights 
w = w*(b-a)/2; 
end

