function tau_d = fiber_delay(distance)
    % fiber_delay calculates the transmission delay through SSMF fiber.
    % Input: 
    %   distance - length of the fiber in meters (L)
    % Output:
    %   tau_d - delay in seconds

    % Constants
    c = 3e8;           % Speed of light in vacuum (m/s)
    n_fiber = 1.468;   % Refractive index of the fiber (typical for silica)
    
    % Velocity of light in the fiber
    v_fiber = c / n_fiber;
    
    % Calculate delay
    tau_d = distance / v_fiber;
end
