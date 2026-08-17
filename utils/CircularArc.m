function [w_top, w_center, w_bottom, theta, s_end] = CircularArc(type, H, H_shell, R,...
    DipAngle, varargin)

%%% This function computes the circular arc slab geometry described used by Johnson et
%%% al., (2017) and Howell & Pappalardo (2019). To run simulations with 
%%% IcySubduction.m using this geometry, set: GeoFlag = 'Johnson' in 
%%% IcySubduction_Parameters.m

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
        [w_top, w_center, w_bottom, theta] = SlabGeometry(s, H, R, DipAngle);
        % Slab = [w_top, w_center, w_bottom];

%%% Compute the cartesian coordinates of the slab column for a given arc length on the
%%% center line.
    case 'location'    
        s = varargin{1};
        [w_top, w_center, w_bottom, theta] = SlabGeometry(s, H, R, DipAngle);
        % Slab = [w_top, w_center, w_bottom];

end

%%% This function computes the difference in depth between the base of the ice shell 
%%% and the bottom surface of the slab at a given arc length s.
    function base = SlabBase(s)

        [~, ~, w_base] = SlabGeometry(s, H, R, DipAngle);
        base = imag(w_base) - H_shell;
    
    end


function [w_top, w_center, w_bottom, theta] = SlabGeometry(s, H, R, DipAngle)

    DipAngle = DipAngle*(pi/180);                           % final slab dip [radians]

%%% Arc length where the slab center line attains the desired dip angle.
    s_bend = DipAngle*R;
    i_r1 = s <= s_bend;
    i_r2 = i_r1 == false;

%%% Dip angle along the slab.
    theta = s/R;
    theta(i_r2) = DipAngle;

%%% Coordinates of the center line for s < R.
    w_center = 1i*R*(1 - exp(1i*s/R)) + 1i*H/2;

%%% Coordinates of the center line for s > R.
    w_center(i_r2) = 1i*R*(1 - exp(1i*DipAngle))...
        + exp(1i*DipAngle)*(s(i_r2) - R*DipAngle) + 1i*H/2;

%%%--------------------------------------------------------------------------------%%%
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