% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
function d = PW_SingleMode_CauchyData(k,l)
% Computes the expansion coefficients of the l-th Fourier mode of the Dirichlet and
% Neumann trace of the plane wave x -> exp(i*k*x_1)
% k = wave number 
% l = integer, Fourier mode number
% This function relies on the Jacobi-Anger expansion of a plane wave into 
% spherical waves. 
% returns the two Fourier coefficients as a column vector 
    
% Bessel function 
    jl =  besselj(l,k);
% Derivative of Bessel function using https://dlmf.nist.gov/10.6, 10.6.2
 if (l == 0)
     jlp = -besselj(1,k);
 else 
     jlp = (besselj(l-1,k) - besselj(l+1,k))/2.0; %J_l'(k)
 end
% First component = Dirichlet mode, second component = Neumann mode    
d = [ (1i^l)* jl; k*(1i^l)*jlp];

end % function d = PW_SingleMode_CauchyData(k,l)
