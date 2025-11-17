function [D] = three_point_centered_varcoeff_D2(z0, zL, n, K, coeff)

%%% function three_point_centered_varcoeff_D2 returns the differentiation 
%%% matrix for computing the second derivative, (K*u_z)z, of a variable u 
%%% over the spatial domain z0 < z < zL. The coefficients, K, vary with z
%%% as a function of u, and so the coefficients of the differentiation
%%% matrix depend on the values of K(z).
%%%
%%%  argument list
%%%
%%%     z0      left value of the spatial independent variable (input)
%%%
%%%     zL      right value of the spatial independent variable (input)
%%%
%%%     n       number of spatial grid points, including the end
%%%             points (input)

%%% Compute the spatial increment.
    dz = (zL - z0)/(n-1);
    
%%% Index vectors for the diagonals.
    e = 2:n-1;
    e_Sub = 2:n;
    e_Sup = 1:n-1;

    switch coeff
        case 'linear'
        %%% Main diagonal.
            M_Diag = (1/2)*[K(1,1) + K(2,1); K(e-1,1) + 2*K(e,1) + K(e+1,1);...
                K(n,1) + K(n-1,1)];

        %%% Sub-diagonal.
            M_SubDiag = -(1/2)*[K(e_Sub-1,1) + K(e_Sub,1); 0];

        %%% Super-diagonal.
            M_SupDiag = -(1/2)*[0; K(e_Sup,1) + K(e_Sup+1,1)];

        case 'harmonic'
        %%% Harmonic mean of the coefficients.
            K_hat = 2*K(1:n-1).*K(2:n)./(K(1:n-1) + K(2:n));

        %%% Main diagonal.
            M_Diag = [K_hat(1); K_hat(e-1,1) + K_hat(e,1);  K_hat(n-1,1)];
    
        %%% Sub-diagonal.
            M_SubDiag = -[K(e_Sub-1,1); 0];
    
        %%% Super-diagonal.
            M_SupDiag = -[0; K(e_Sup+1,1)];
    end
    
%%% For the differentiation operator.
    M = spdiags([M_SubDiag M_Diag M_SupDiag], -1:1, n, n);
    D = -M/dz^2;
    
%%% Summation-by-parts form of the operator from Erickson and Dunham (2014). It's off by a
%%% factor of 1/dz. But the equations here look equivalent to those in Erickson2014, so
%%% perhaps there's a typo in Erickson's equations.
    % H = (dz/2)*spdiags([1; 2*ones(n-2,1); 1], 0, n, n);
    % B = spdiags([-K(1,1); zeros(n-2,1); K(n,1)], 0, n, n);    
    % S = (1/dz)*spdiags([[zeros(n-3,1);1/2;0;0] [zeros(n-2,1);-2;0]...
    %     [-3/2;ones(n-2,1);3/2] [0;2;zeros(n-2,1)] [0;0;-1/2;zeros(n-3,1)]],...
    %     -2:2, n, n);    
    % D = H^(-1)*(-M + B*S);
    