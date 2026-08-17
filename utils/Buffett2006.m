function [w_top, w_center, w_bottom, theta, K, dKds, s_end] = Buffett2006(type, H, ...
    H_shell, R_min, varargin)

%%% This function computes the slab geometry described in paragraph 10 of Buffett 
%%% (2006). The geometry is apparently based on observations of Earth subduction 
%%% zones. Although no references are cited in that paragraph of Buffett (2006) , they 
%%% can be found in Buffett & Becker (2012).

%%% Buffett assumes that the plate begins bending at s = 0, begins unbending at
%%% s = R_min, and becomes straight at s = 2*R_min. Here, s is the arc length along
%%% the slab's mid-surface at z = 0. Buffett refers to the mid-surface as the 
%%% centerline.

%%% Buffett, B. A. (2006), Plate force due to bending at subduction zones, 
%%% J. Geophys. Res., 111, B09405, doi:10.1029/2006JB004295.

%%% Buffett, B. A., and T. W. Becker (2012), Bending stress and dissipation in 
%%% subducted lithosphere, J. Geophys. Res., 117, B05413, doi:10.1029/2012JB009205.

%%%-------------------------------------------------------------------------------%%%

switch type

%%% Output the entire slab geometry by computing the value of the centerline arc 
%%% length that corresponds to where the bottom surface of the slab intersects the
%%% base of the ice shell. This is used for plotting the geometry and also computing 
%%% the run time for subduction simulations.
    case 'geometry'

%%% The total arc length is determined iteratively.
        s_end = fzero(@SlabBase, H_shell);

%%% Now compute the complete slab geometry.
        s = linspace(0, s_end, 100)';
        [w_top, w_center, w_bottom, theta, K, dKds] = SlabGeometry(s, H, R_min);
        % Slab = [w_top, w_center, w_bottom];

%%% Compute the cartesian coordinates of the slab column for a given arc length on 
%%% the center line.
    case 'location'    
        s = varargin{1};
        [w_top, w_center, w_bottom, theta, K, dKds] = SlabGeometry(s, H, R_min);
        % Slab = [w_top, w_center, w_bottom];

end

%%% This function computes the difference in depth between the base of the ice shell 
%%% and the bottom surface of the slab at a given arc length s.
    function base = SlabBase(s)

        [~, ~, w_base] = SlabGeometry(s, H, R_min);
        base = imag(w_base) - H_shell;
    
    end

function [w_top, w_center, w_bottom, theta, K, dKds] = SlabGeometry(s, H, R_min)

%%% Allocate arrays.
    dKds = nan(size(s));            % gradient of curvature along center line.
    K = nan(size(s));               % curvature of center line.
    theta = nan(size(s));           % local dip angle
    w_center = nan(size(s));        % complex Cartesian coordinates of center line.

%%% Use the analytic expressions to compute the complex coordintes (w = x + iy) of 
%%% the slab center line and the slab dip angle. Note that the dip angle is constant 
%%% along the z-axis of the slab.
    i_r1 = s <= R_min;
    i_r2 = s > R_min & s <= 2*R_min;
    i_r3 = s > 2*R_min;

%%% Geometry for s < R_min.
    dKds(i_r1) = 1/R_min^2;
    K(i_r1) = s(i_r1)/R_min^2;
    theta(i_r1) = s(i_r1).^2/(2*R_min^2);
    w_center(i_r1) = -(-1)^(3/4)*sqrt(pi/2)*R_min*erfi((-1)^(1/4)*s(i_r1)/(sqrt(2)*R_min));
    

%%% Geometry for R_min < s < 2*R_min.
    a = 1/(2*R_min^2);
    b = 2/R_min;
    dKds(i_r2) = -1/R_min^2;
    K(i_r2) = b - 2*a*s(i_r2);
    
    theta(i_r2) = b*s(i_r2) - a*s(i_r2).^2 - 1;
    
    w_center(i_r2) = -(-1)^(3/4)*sqrt(pi/2)*R_min*erfi((-1)^(1/4)/sqrt(2))...
        -(-1)^(1/4)*(pi/(4*a))^(1/2)*exp(1i*b^2/(4*a) - 1i)...
        *(erfi((-1)^(3/4)*(2*a*s(i_r2) - b)/(2*sqrt(a)))...
        - erfi((-1)^(3/4)*(2*a*R_min - b)/(2*sqrt(a))));

%%% Coordinates and dip angle for s > 2*R_min.
    dKds(i_r3) = 0;
    K(i_r3) = 0;
    theta(i_r3) = 1;
    w_center(i_r3) = -(-1)^(3/4)*sqrt(pi/2)*R_min*erfi((-1)^(1/4)/sqrt(2))...
        -(-1)^(1/4)*(pi/(4*a))^(1/2)*exp(1i*b^2/(4*a) - 1i)...
        *(erfi((-1)^(3/4)*(2*a*2*R_min - b)/(2*sqrt(a)))...
        - erfi((-1)^(3/4)*(2*a*R_min - b)/(2*sqrt(a))))...
        + exp(1i).*(s(i_r3) - 2*R_min);

%%% Shift the center line to half the slab thickness.
    w_center = w_center + 1i*H/2;

%%%-------------------------------------------------------------------------------%%%
%%% Find the coordinates of the top and bottom surfaces of the slab by projecting a 
%%% distance H/2 along the angle normal to the center line.

%%% Normal angle (psi) along the center line.
    psi = pi/2 - theta;
    x_center = real(w_center);
    y_center = imag(w_center);

%%% Coordinates of the top and bottom surfaces.
    dx = (H/2)*cos(psi);
    dy = (H/2)*sin(psi);

    x_top = x_center + dx;
    x_bottom = x_center - dx;

    y_top = y_center - dy;
    y_bottom = y_center + dy;
    
    w_top = x_top + 1i*y_top;
    w_bottom = x_bottom + 1i*y_bottom;
end   

end