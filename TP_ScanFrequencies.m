% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
function [l2n,h1n] = TP_ScanFrequencies(k_range,ni,Lmax)
% Compute and plot the norms of the solution for the Helmholtz transmission problem for a
% range of frequencies
    
    % Default number of Fourier modes
    if (nargin < 3), Lmax = 50; end 
    
    l2n = []; h1n = []; 
    for k=k_range
        [l2k,h1k] = PW_Scattering_SolNorms(k,ni,Lmax);
        fprintf('TP_ScanFrequencies: k = %d, L2 norm = %d, h1-norm = %d\n',k,l2k,h1k);
        l2n = [l2n;l2k];
        h1n = [h1n,h1k];
    end 
    
    figure('name','Norms');
    plot(k_range,l2n,'b-',k_range,h1n,'r-');
    xlabel('wave number k','fontsize',14);
    ylabel('solution norms in ball B_2','fontsize',14);
    title(sprintf('{Helmholtz T.P., incident p.w., n_i = %f}',ni));
    legend('L2 norm','H1 norm','location','best');
end
