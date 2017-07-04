    %THIS CODE REPLICATES THE SOME OF THE RESULTS FROM HE-LITTERMAN 1999
clear;
clc;
%INPUT CORRELATION, VOL, AND E[R] MATRICES
rho = [ 1.000 0.488 0.478 0.515 0.439 0.512 0.491;
      0.488 1.000 0.664 0.655 0.310 0.608 0.779;
      0.478 0.664 1.000 0.861 0.355 0.783 0.668;
      0.515 0.655 0.861 1.000 0.354 0.777 0.653;
      0.439 0.310 0.355 0.354 1.000 0.405 0.306;
      0.512 0.608 0.783 0.777 0.405 1.000 0.652;
      0.491 0.779 0.668 0.653 0.306 0.652 1.000];

weq = [0.016,0.022,0.052,0.055,0.116,0.124,0.615];
sigma = [0.160 0.203 0.248 0.271 0.210 0.200 0.187];
assets={'Australia';'Canada   ';'France   ';'Germany  ';'Japan    ';'UK       ';'USA      '};
cov = (sigma' * sigma) .* rho;

%WE WANT TO OPTIMIZE W'U-.5*DELTA*W'*SIGMA*W 
delta = 2.5; %PAGE 10 OF PAPER
tau = .05; %PAGE 11 OF PAPER
EqRiskPremia = delta * cov *  weq'; %PAGE 3 OF PAPER
    
%GERMAN WILL OUTPERFORM EUROPE BY 5%
P = [0 0 -.295 1 0 -.705 0]; %MARKET CAP WEIGHT FRANCE AND UK (E.G. FRANCE = .052 / (.052 + .124);
Q = .05;
ts = tau * cov;
omega = P * ts * P' .* eye(rows(Q));
mbar_inv = inv(inv(ts) + P' * inv(omega) * P );
sigma_bar = cov + mbar_inv;
mu_bar = mbar_inv * (inv(ts) * EqRiskPremia + P' * inv(omega) * Q );
weight = inv(delta .* sigma_bar) * mu_bar;

Table4Col4 = weight - (weq') / (1 + tau);
%PRINT OUT TABLE THAT MATCHES TABLE 4 OF THE PAPER
Table4 = [{'Assets' 'P' 'Mu_bar' 'Opt. Weights' 'Deviation'};
            assets num2cell([P' mu_bar weight ...
            round(Table4Col4 * 1000)./1000] * 100)]

%CANADA OUTPERFORMS US BY 3% A YEAR;
P = [0 0 -.295 1 0 -.705 0; 0 1 0 0 0 0 -1];
Q = [.05; .03];
omega = P * ts * P' .* eye(rows(Q));
mbar_inv = inv(inv(ts) + P' * inv(omega) * P );
sigma_bar = cov + mbar_inv;
mu_bar = mbar_inv * (inv(ts) * EqRiskPremia  + P' * inv(omega) * Q );
weight = inv(delta * sigma_bar) * mu_bar;

Table5Col5 = weight - (weq') / (1 + tau);
%PRINT OUT TABLE THAT MATCHES TABLE 4 OF THE PAPER
Table5 = [{'Assets' 'P-1' 'P-2' 'Mu_bar' 'Opt. Weights' 'Deviation'};
            assets num2cell([P' mu_bar weight ...
            round(Table5Col5 * 1000)./1000] * 100)]
%PRECISION RATIO W/TAU AND WEIGHT ON EACH VIEW        
Confidence_ratio = diag(omega) ./ tau;
Lambda = pinv(P)' * (weight * (1 + tau) - weq');

