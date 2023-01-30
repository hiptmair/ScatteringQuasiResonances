% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
function A = OpMat_TP_SolOp(k,l,n_i,n_o)
% Matrix representation of the solution operator for Fourier mode l and with respect to
% Bessel and Hankel function basis for the inner and outer subdomains, see notes
% ScatteringAtSphere for details. 
% k = wave number 
% l = number of Fourier mode, non-negative 
% n_i, n_o = inner and outer refractive index 
% The function returns a 2x2 matrix describing the mapping of the two l-th Fourier 
% coefficients of Cauchy (Dirichlet and Neumann) jump data to the corresponding expansion
% coefficients for the inner and outer solution
% Note: mapping matrix depends only on the modulus of the mode number 
    
    % Effective wave numbers 
    k_i = sqrt(n_i)*k;
    k_o = sqrt(n_o)*k;

    % Bessel function  J_l(k_i)
    jl_i = besselj(l,k_i); 
    % Hankel function of the first kind H^1_l(k_o)
    hl_o = besselh(l,k_o); 
    % Derivatives of Bessel and Hankel functions
    if (l == 0)
        % See https://dlmf.nist.gov/10.6, 10.6.3
         jlp_i = -besselj(1,k_i);
         hlp_o = - besselh(l+1,k_o);
    else 
        % We use the well-known formulas expressing the the derivatives of 
        % those functions as difference of the same function with orders
        % raised/lowered by one, see https://dlmf.nist.gov/10.6, 10.6.2
       jlp_i = (besselj(l-1,k_i) - besselj(l+1,k_i))/2.0; %J_l'(k_i)
       hlp_o = (besselh(l-1,k_o) - besselh(l+1,k_o))/2.0; % H^1_l'(k_o)
    end 
    % Generate mapping matrix 
    A = [hl_o      , -jl_i; ...
         k_o*hlp_o , -k_i*jlp_i];
end % function A = OpMat_TP_SolOp(k,l,n_i,n_o)


