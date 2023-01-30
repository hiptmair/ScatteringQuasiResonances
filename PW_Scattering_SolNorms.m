% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
function [l2norm,h1norm] = PW_Scattering_SolNorms(k,ni,Lmax)
% Computes the solution of a Helmholtz transmission problem with wave number k, refective
% index 1 outside the unit disk, refractive index ni inside the unit disk, and excitation
% by a plane wave x -> exp(ik*x_1) by means of Fourier spectral analysi using the lowest 
% Lmax Fourier modes. The function returns the L2 norm and H1 norm of the solution in a
% disk of radius 2
    
    fprintf('PW_Scattering_SolNorms: k=%d, ni = %d, Lmax = %d\n',k,ni,Lmax); 
    
    % Compute Fourier coefficients of Cauchy data for incident plance wave (Jacobi-Anger
    % expansion)
    coeff_uinc = zeros(2,2*Lmax+1);
    for l=(-Lmax:Lmax)
        idx = l+Lmax+1;
        coeff_uinc(:,idx) = PW_SingleMode_CauchyData(k,l);
    end

    % Solve for every single mode. Store expansion coefficients of the solution in a 2xL
    % matrix, top row = exterior expansion coefficients, bottom row = interior expansion coefficients
    coeff_sol = zeros(2,2*Lmax+1);
    for l=(-Lmax:Lmax)
        idx = l+Lmax+1;
        A = OpMat_TP_SolOp(k,abs(l),ni,1.0);
        coeff_sol(:,idx) = A\coeff_uinc(:,idx);
    end

    % Compute norms of the solution
    [l2norm,h1norm] = TP_Sol_Norms(k,ni,1.0,coeff_sol(2,:),coeff_sol(1,:));
end
