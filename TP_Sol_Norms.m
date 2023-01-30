% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
function [l2norm,h1norm] = TP_Sol_Norms(k,n_i,n_o,ci_vec,co_vec)
% Computes the L2 and H1 norm of a solution of the Helmholtz transmission problem for a
% unit-disk obstacle over the disk with radius = 2
% k = wave number 
% n_i,n_o = interior and exterior refractive indices 
% ci_vec = vector of expansion coefficients for interior solution
% co_vec = vector of expansion coefficients for exterior solution 
% These two vectors need not have the same length
    
    % Coefficient vectors should be column vectors 
    ci_vec = reshape(ci_vec,length(ci_vec),1);
    co_vec = reshape(co_vec,length(co_vec),1);

    % Obtain nodes and weight for Gaussian quadrature on [0,1]
    n_gauss_pts = 200;
    [x,w] = gaussquad(0.0,1.0,n_gauss_pts);
    
    % Effective wave numbers 
    k_i = sqrt(n_i)*k;
    k_o = sqrt(n_o)*k;

    % Interior solution, coefficient vector ci_vec   
    % Find index range for Fourier coefficients
    Ni = length(ci_vec); mi = floor(Ni/2);
    if (mod(Ni,2) == 1)
        idxi_min = -mi;
        idxi_max = mi;
    else
        idxi_min = -mi+1;
        idxi_max = mi;
    end
    % Integral of |J_l(k_i*r)|^2*r over [0,1], 
    % computed with Gauss quadrature for l=0,...mi
    L2IJl = zeros(mi+1,1);
    for l=(0:mi) 
        L2IJl(l+1) = dot(w,abs(besselj(l,k_i*x)).^2.*x);
    end
    L2IJl = [L2IJl(-idxi_min+1:-1:2);L2IJl];
    % Integral of |k_i*J_l'(r)(k_i*r)|^2*r + l/r|J_l(k_i*r)|^2
    % Special case: l=0
    H1IJl =  zeros(mi+1,1);
    H1IJl(1) = dot(w,abs(k_i*besselj(1,k_i*x)).^2.*x);
    for  l=(1:mi) 
        Jlp = @(r) ((besselj(l-1,r) - besselj(l+1,r))/2.0);
        H1IJl(l+1) = dot(w,abs(k_i*Jlp(k_i*x)).^2.*x + l*l*(abs(besselj(l,k_i*x)).^2)./x);
    end
    H1IJl = [H1IJl(-idxi_min+1:-1:2);H1IJl];
    % Compute norms over unit disk
    l2norm_i_sq = 2*pi*dot(abs(ci_vec).^2,L2IJl);
    h1snorm_i_sq = 2*pi*dot(abs(ci_vec).^2,H1IJl);
    
    % Exterior solution defined by coefficient vector co_vec
    % Obtain nodes and weight for Gaussian quadrature on [1,2]
    n_gauss_pts = 200;
    [x,w] = gaussquad(1.0,2.0,n_gauss_pts);
    % Find index range for Fourier coefficients
    No = length(co_vec); mo = floor(No/2);
    if (mod(No,2) == 1)
        idxo_min = -mo;
        idxo_max = mo;
    else
        idxo_min = -mo+1;
        idxo_max = mo;
    end
    % Integral of |H^1_l(k_o*r)|^2*r over [1,2], 
    % computed with Gauss quadrature for l=0,...mi
    L2IHl = zeros(mo+1,1);
    for l=(0:mo) 
        L2IHl(l+1) = dot(w,abs(besselh(l,k_o*x)).^2.*x);
    end
    L2IHl = [L2IHl(-idxo_min+1:-1:2);L2IHl];
    % Integral of |k_i*H^1_l'(r)(k_o*r)|^2*r + l/r|H^1_l(k_o*r)|^2
    % Special case: l=0
    H1IHl =  zeros(mo+1,1);
    H1IHl(1) = dot(w,abs(k_o*besselh(1,k_o*x)).^2.*x);
    for  l=(1:mo) 
        Hlp = @(r) ((besselj(l-1,r) - besselj(l+1,r))/2.0);
        H1IJl(l+1) = dot(w,abs(k_o*Hlp(k_o*x)).^2.*x + l*(abs(besselh(l,k_o*x)).^2)./x);
    end
    H1IHl = [H1IHl(-idxo_min+1:-1:2);H1IHl];
    % Compute norms over the annulus
    l2norm_o_sq = 2*pi*dot(abs(co_vec).^2,L2IHl);
    h1snorm_o_sq = 2*pi*dot(abs(co_vec).^2,H1IHl);
    % Compute final norms
    l2norm = sqrt(l2norm_i_sq + l2norm_o_sq);
    h1norm = sqrt(l2norm_i_sq + l2norm_o_sq + h1snorm_i_sq + h1snorm_o_sq);
end
