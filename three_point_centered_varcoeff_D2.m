function [D] = three_point_centered_varcoeff_D2(z0,zL,n,K)

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
    
%    D = zeros(n,n);
%    iD = (1:N);
%    iSup = (1:N-1);
%    iSupSup = (1:N-2);
    
%%% Interior nodes.
    e = 2:n-1;    
    M_Diag = (1/2)*[K(1,1) + K(2,1); K(e-1,1) + 2*K(e,1) + K(e+1,1);...
        K(n,1) + K(n-1,1)];
    
    eSub = 2:n;
    M_SubDiag = -(1/2)*[K(eSub-1,1) + K(eSub,1); 0];
    
    eSup = 1:n-1;
    M_SupDiag = -(1/2)*[0; K(eSup,1) + K(eSup+1,1)];
    
    M = spdiags([M_SubDiag M_Diag M_SupDiag], -1:1, n, n);
    
    H = (dz/2)*spdiags([1; 2*ones(n-2,1); 1], 0, n, n);
    
    B = spdiags([-K(1,1); zeros(n-2,1); K(n,1)], 0, n, n);
    
    S = (1/dz)*spdiags([[zeros(n-3,1);1/2;0;0] [zeros(n-2,1);-2;0]...
        [-3/2;ones(n-2,1);3/2] [0;2;zeros(n-2,1)] [0;0;-1/2;zeros(n-3,1)]],...
        -2:2, n, n);
    
    D = H^(-1)*(-M + B*S);