function y_top = Buffett2006_Top(s_end, H, H_shell, R_min, x_Slab)

%%% For a given location in the slab (x_0, y_0), find the location of the top of the 
%%% slab that is directly above. I.E. find y_top(x_0).

%%%-------------------------------------------------------------------------------%%%
%%% Find the zero of the function for the x-coordinate of the top surface.
    y_top = nan(size(x_Slab));

    for i = 1:numel(x_Slab)
        x_0 = x_Slab(i);
        fun = @f_top;
        x_int = [0 2*s_end];
        % try
            s_0 = fzero(fun, x_int);
        % catch
        %     i
        % end

%%% Evaluate the geometry to get the y-coordinate of the top surface.
        type = 'location';
        [w_top, ~, ~, ~, ~, ~] = Buffett2006(type, H, H_shell, R_min, s_0);
        y_top(i) = imag(w_top);
    end

function zero_top = f_top(s)
    
%%% Use the analytic expressions to compute the complex coordintes (w = x + iy) of 
%%% the slab center line and the slab dip angle. This code contains the same 
%%% equations as Buffet2006.m.
    a = 1/(2*R_min^2);
    b = 2/R_min;

    if s <= R_min
%%% Geometry for s < R_min.
        theta = s.^2/(2*R_min^2);

        w_center = -(-1)^(3/4)*sqrt(pi/2)*R_min*erfi((-1)^(1/4)*s/(sqrt(2)*R_min));
    
    elseif s > R_min && s <= 2*R_min
%%% Geometry for R_min < s < 2*R_min.
        theta = b*s - a*s.^2 - 1;

        w_center = -(-1)^(3/4)*sqrt(pi/2)*R_min*erfi((-1)^(1/4)/sqrt(2))...
            -(-1)^(1/4)*(pi/(4*a))^(1/2)*exp(1i*b^2/(4*a) - 1i)...
            *(erfi((-1)^(3/4)*(2*a*s - b)/(2*sqrt(a)))...
            - erfi((-1)^(3/4)*(2*a*R_min - b)/(2*sqrt(a))));

    elseif s > 2*R_min
%%% Coordinates and dip angle for s > 2*R_min.
        theta = 1;
        w_center = -(-1)^(3/4)*sqrt(pi/2)*R_min*erfi((-1)^(1/4)/sqrt(2))...
            -(-1)^(1/4)*(pi/(4*a))^(1/2)*exp(1i*b^2/(4*a) - 1i)...
            *(erfi((-1)^(3/4)*(2*a*2*R_min - b)/(2*sqrt(a)))...
            - erfi((-1)^(3/4)*(2*a*R_min - b)/(2*sqrt(a))))...
            + exp(1i).*(s - 2*R_min);
    end

%%% Shift the center line to half the slab thickness.
    w_center = w_center + 1i*H/2;

%%%-------------------------------------------------------------------------------%%%
%%% Find the coordinates of the top and bottom surfaces of the slab by projecting a 
%%% distance H/2 along the angle normal to the center line.

%%% Normal angle (psi) along the center line.
    psi = pi/2 - theta;
    x_center = real(w_center);

%%% x-coordinate of the top surface.
    dx = (H/2)*cos(psi);
    x_top_s = x_center + dx;

%%%-------------------------------------------------------------------------------%%%
%%% The desired solution is the zero of x_top - x_top_s.
    zero_top = x_0 - x_top_s;
end

end