%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the example code for the Recursive Trigonometry (RT) technique,
% which demonstrates the use of recursive methods to compute trigonometric
% functions (e.g., cosine) over a full period with reduced computational
% complexity. The code also compares the RT technique with the Radix-2 and
% Radix-4 Coordinate Rotation Digital Computer (CORDIC) algorithms to
% evaluate their computational accuracy and efficiency.
%
% The main objectives of this code are:
% 1. Implement the Recursive Trigonometry (RT) technique to calculate cosine
%    values for a full cycle (0 to 2*pi).
% 2. Perform iterative computations for sine and cosine using the Radix-2
%    and Radix-4 CORDIC algorithms.
% 3. Analyze the performance of these methods by calculating and plotting
%    the error deviations between the computed and reference values.
% 4. Compute Mean Squared Error (MSE) and Root Mean Squared Error (RMSE) to
%    quantify the accuracy of the techniques.
%
% This example serves as a foundational reference for applications that
% require efficient trigonometric computations, such as digital signal
% processing, communication systems, and embedded systems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

% Define the step size Seta for angle increments in radians
Seta = 0.0062831853;

% Calculate the number of points for specific angular ranges based on Seta
Num_1 = floor(pi / (2 * Seta));        % Number of points for 90 degrees (pi/2)
Num_2 = floor(pi / Seta);             % Number of points for 180 degrees (pi)
Num_3 = floor(3 * pi / (2 * Seta));   % Number of points for 270 degrees (3*pi/2)
Num_4 = floor(2 * pi / Seta);         % Number of points for 360 degrees (2*pi)
Num = Num_2;                          % Number of iterations up to 180 degrees (pi)

%% RT Technique
test = vpa(cos(Seta), 7);             % Calculate cos(Seta) with 7-digit precision
C0 = test;                            % Assign cos(Seta) to C0

C = 2 * C0;                           % Precompute 2 * cos(Seta)
A2 = C * C0 - 1;                      % Calculate cos(2*Seta) using a recursive formula
A1 = C0;                              % Initialize A1 as cos(Seta)

new_out(1) = C0;                      % Store cos(Seta) as the first value in the result array
new_out(2) = A2;                      % Store cos(2*Seta) as the second value in the result array

coss(1) = cos(Seta);                  % Reference (ground truth) for cos(Seta)
coss(2) = cos(2 * Seta);              % Reference (ground truth) for cos(2*Seta)

% Loop to compute cos(i*Seta) for one full period
for i = 3:Num                         % Calculate cosine values for one period
    if i <= Num_1
        % First quadrant (0 to pi/2)
        A3 = C * A2 - A1;             % Compute cos(i*Seta) recursively
        A1 = A2;                      % Update A1 for the next iteration
        A2 = A3;                      % Update A2 for the next iteration
        new_out(i) = A3;              % Store the computed value
        if i == Num_1
            m = 0;                    % Reset symmetry index for the second quadrant
        end
        coss(i) = cos(i * Seta);      % Store ground truth cosine value
    elseif (i > Num_1) & (i <= Num_2)
        % Second quadrant (pi/2 to pi)
        new_out(i) = -new_out(Num_1 - m);   % Use symmetry property: cos(pi-x) = -cos(x)
        coss(i) = -cos((Num_1 - m) * Seta); % Ground truth using symmetry
        m = m + 1; 
        if i == Num_2
            m = 0;                    % Reset symmetry index for the third quadrant
        end
    elseif (i > Num_2) & (i <= Num_3)
        % Third quadrant (pi to 3*pi/2)
        new_out(i) = -new_out(m + 1);       % Use symmetry property: cos(pi+x) = -cos(x)
        coss(i) = -cos((m + 1) * Seta);     % Ground truth using symmetry
        m = m + 1;
        if i == Num_3
            m = 0;                    % Reset symmetry index for the fourth quadrant
        end 
    elseif (i > Num_3) & (i <= Num_4)
        % Fourth quadrant (3*pi/2 to 2*pi)
        new_out(i) = new_out(Num_1 - m);    % Use symmetry property: cos(2*pi-x) = cos(x)
        coss(i) = cos((Num_1 - m) * Seta);  % Ground truth using symmetry
        m = m + 1; 
        if i == Num_4
            m = 0;    
        end       
    end
    % Calculate the absolute error between the computed and ground truth cosine values
    NE(i) = abs(coss(i) - new_out(i));  
end

% Compute Mean Squared Error (MSE) and Root Mean Squared Error (RMSE)
D = NE.^2;
MSE = sum(D(:)) / numel(coss);
RMSE = sqrt(MSE);

% Define angular range for plotting
g = 0:Seta:(Num * Seta - Seta);

% Plot the error for the RT Technique
subplot(3,1,3)
plot(g, NE, 'r-');                    % Plot the error curve
xticks([0 pi/2 pi]);                  % Set x-axis tick positions
xticklabels({'0', '90', '180'});      % Label x-axis ticks in degrees
xlabel('Angle (Degree)', 'FontSize', 15);   % Label for x-axis
ylabel('Error (Deviation)', 'FontSize', 15); % Label for y-axis
set(gca, 'YLim', [0 1.6*10^-11], 'FontSize', 15); % Set y-axis limits
text(-0.15, 0.5, '(c)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15);

% Output the maximum error for the RT technique
max_NE = max(NE);
fprintf('Max value for NE: %.2e\n', max_NE);

%% radix-2 CORDIC
clc;

% Initialize scaling factor and angular step size
k = 0.607253;                   % Scaling factor for CORDIC
angle = 0;                      % Initialize angle
Seta = 0.0062831853;            % Define angular step size

Num = floor(pi / Seta);         % Number of iterations up to pi

%%%%%%%%%%%%%%%%%%%%%%    CORDIC  %%%%%%%%%%%%%%%%%%%%%%%%%
angle = 0;
angle1 = 0;

for i = 1:Num
    angle = angle + Seta; 
    if angle >= 2*pi     
        angle = 0;               % Reset angle after a full period
    end
    
    % Map angle to the first quadrant
    if angle >= pi*3/2
        z0 = angle - pi*3/2;
    elseif angle >= pi
        z0 = angle - pi;
    elseif angle >= pi/2
        z0 = angle - pi/2;
    else
        z0 = angle;
    end  
    
    % Initialize variables for the iterative CORDIC calculation
    x0 = k;
    y0 = 0;
    z = z0;
    a = 0;
  
    % Perform 17 iterations to compute sine and cosine values
    while a < 17
        if z0 > 0
            d = 1;               % Positive rotation direction
        else
            d = -1;              % Negative rotation direction
        end  
        
        % Update x, y, and z using the CORDIC formulas
        x1 = x0 - d * y0 * 2^(-a);
        y1 = y0 + d * x0 * 2^(-a);
        z1 = z0 - d * atan(2^(-a));
        
        % Update variables for the next iteration
        x0 = x1;
        y0 = y1;
        z0 = z1;
        a = a + 1;
    end
    x = x1;                      % Final x-coordinate (cosine)
    y = y1;                      % Final y-coordinate (sine)

    % Map computed values back to the appropriate quadrant
    if angle >= pi*3/2
        sin_out(i) = -x;
        cos_out(i) = y;
    elseif angle >= pi
        sin_out(i) = -y;
        cos_out(i) = -x;
    elseif angle >= pi/2
        sin_out(i) = x;
        cos_out(i) = -y;
    else
        sin_out(i) = y;
        cos_out(i) = x;
    end
    rr = cos_out(i);             % Assign computed cosine value to rr
    
    % Calculate the absolute error between the computed and reference cosine values
    Right = cos(angle);    
    CD(i) = abs(Right - rr);  
end

% Compute MSE and RMSE for the CORDIC method
D1 = CD.^2;
MSE1 = sum(D1(:)) / numel(Right);
RMSE1 = sqrt(MSE1);

%% radix-4 CORDIC
clc;

% Number of iterations (die) and initialize arrays
die = 16;                          % Number of iterations for Radix-4
xs = zeros(die+1, 1);              % Array to store x-coordinates
ys = zeros(die+1, 1);              % Array to store y-coordinates
zs = zeros(die+1, 1);              % Array to store angles

Setas = 0.0062831853;              % Step size in radians
Nums = floor(pi / Setas);          % Number of steps for half a period
angle5 = 0;                        % Initialize angle accumulator

% Loop over the number of steps (Nums) to compute cosine values
for n = 1:Nums
    angle5 = angle5 + Setas;       % Increment angle by step size
    if angle5 >= 2*pi
        angle5 = 0;                % Reset angle after a full period
    end
    
    % Map the angle to the first quadrant
    if angle5 >= pi*3/2
        z2 = angle5 - pi*3/2;
    elseif angle5 >= pi
        z2 = angle5 - pi;
    elseif angle5 >= pi/2
        z2 = angle5 - pi/2;
    else
        z2 = angle5;
    end  
    
    % Initialize CORDIC variables for the first iteration
    x2 = 1;                        % Initial x-coordinate
    y2 = 0;                        % Initial y-coordinate
    w = z2;                        % Initialize angle

    % Determine direction (ds) for Radix-4 based on initial angle
    if w >= (5/8)
        ds = 2;
    elseif (w >= 3/8) & (w < 5/8)
        ds = 1;
    elseif (w >= -1/2) & (w < 3/8)
        ds = 0;
    elseif (w >= -7/8) & (w < -1/2)
        ds = -1;
    else
        ds = -2;
    end

    % Perform the first iteration
    k0 = (1 + ds^2)^(-1/2);         % Scaling factor for the first iteration
    xs(1) = x2 - ds*y2;             % Update x-coordinate
    ys(1) = y2 + ds*x2;             % Update y-coordinate
    zs(1) = z2 - atan(ds);          % Update angle

    % Perform subsequent iterations
    for is = 1:die
        w = 4^(is) * zs(is);        % Scale the angle for the current iteration

        % Determine direction (ds) based on scaled angle
        if w >= (3/2)
            ds = 2;
        elseif (w >= 1/2) & (w < 3/2)
            ds = 1;
        elseif (w >= -1/2) & (w < 1/2)
            ds = 0;
        elseif (w >= -3/2) & (w < -1/2)
            ds = -1;
        else
            ds = -2;
        end

        % Update x, y, and z values using Radix-4 formulas
        xs(is+1) = xs(is) - ds*ys(is)*(4^(-is)); % Update x-coordinate
        ys(is+1) = ys(is) + ds*xs(is)*(4^(-is)); % Update y-coordinate
        zs(is+1) = zs(is) - atan(ds*4^(-is));    % Update angle

        % Calculate scaling factor for the current iteration
        k(is) = (1 + (ds^2) * 4^(-2*is))^(-1/2);
    end

    % Compute cumulative scaling factor
    Kv = k0;
    for j = 1:is
        Kv = Kv * k(j);
    end

    % Apply scaling factor to final x and y values
    a = xs(17) * Kv;
    b = ys(17) * Kv;

    % Assign sine and cosine values based on the quadrant
    if angle5 >= pi*3/2
        sinout(n) = -a;
        cosout(n) = b;
    elseif angle5 >= pi
        sinout(n) = -b;
        cosout(n) = -a;
    elseif angle5 >= pi/2
        sinout(n) = a;
        cosout(n) = -b;
    else
        sinout(n) = b;
        cosout(n) = a;
    end

    % Compute the absolute error between computed and true cosine values
    ss(n) = single(cosout(n));      % Convert computed cosine to single precision
    rights = cos(angle5);           % True cosine value
    FV(n) = abs(rights - ss(n));    % Absolute error
end

% Compute MSE and RMSE for Radix-4 CORDIC
D2 = FV.^2;
MSE2 = sum(D2(:)) / numel(rights);
RMSE2 = sqrt(MSE2);

%% Plot
% Define the angular range for plotting
g = 0:0.0062831853:(pi-0.0062831853);

% Plot Radix-2 CORDIC error
subplot(3,1,1)
plot(g, CD, 'b-');                 % Plot Radix-2 error in blue
ylabel('Error (Deviation)', 'FontSize', 15); % Label for y-axis
set(gca, 'YLim', [0 1.6*10^-5], 'FontSize', 15); % Set y-axis limits
text(-0.15, 0.5, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15);
set(gca, 'XTick', [], 'FontSize', 15); % Remove x-axis ticks

% Calculate and display maximum error for Radix-2 CORDIC
max_CD = max(CD);
fprintf('Max value for CD: %.2e\n', max_CD);

% Plot Radix-4 CORDIC error
subplot(3,1,2)
plot(g, FV, 'm-');                 % Plot Radix-4 error in magenta
ylabel('Error (Deviation)', 'FontSize', 15); % Label for y-axis
set(gca, 'YLim', [0 3.5*10^-8], 'FontSize', 15); % Set y-axis limits
text(-0.15, 0.5, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15);
set(gca, 'XTick', [], 'FontSize', 15); % Remove x-axis ticks

% Calculate and display maximum error for Radix-4 CORDIC
max_FV = max(FV);
fprintf('Max value for FV: %.2e\n', max_FV);

