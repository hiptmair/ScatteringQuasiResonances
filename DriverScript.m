% Supplementary material for the article
% "Frequency-Explicit Shape Uncertainty Quantification for Acoustic Scattering
% by R. Hiptmair, Ch. Schwab, and E. A. Spence
%
% Driver script for numerical experiments concerning scattering at a disk
% Data are saved in variables for further exploration
k_range = 2:0.001:13.5;
[l2n_three,h1n_three] = TP_ScanFrequencies(k_range,3.0);
[l2n_third,h1n_third] = TP_ScanFrequencies(k_range,1/3.0);
