function[angular_size] = Pix2Angle(pixel_size, monitor_distance, direction)
%
% This function provides a conversion from linear size of a line segment
% (in pixels) to the angular size subtended by the line (in degrees).
% 
% Direction is a string. It must be 'h' or 'v'. This determines whether
% the line segment is oriented vertically or horizontally.
%
% Monitor distance must be given in cm
%
% THIS CODE IS SPECIFIC TO THE OLD MONITOR IN Bill's RIG BY VIRTUE OF
% HARD CODING THE EFFECTIVE PICTURE SIZE AND NUMBER OF PIXELS
%
% Color Graphic Display GDM-20D11 Manual lists: (Model number has been double
% checked!)
% "Effective picture size" = 384 x 290 mm
% "Display picture size" = 373 x 280 mm.	<---**** I will use this.
%

% Hard codes

SCRXPIX = 1280.0;       % # of horizontal pixels
SCRYPIX = 1024.0;       % # of vertical pixels
HMONCM = 37.3;          % Effective picture horz size 
VMONCM = 28.0;          % Effective picture vert size 


% Step 1: Find the number of pixels per cm on screen

if (direction == 'h'),
        pixels_per_cm = (SCRXPIX / HMONCM);
elseif (direction == 'v'),
        pixels_per_cm = (SCRYPIX / VMONCM);
else
        disp('Error in direction variable in Pix2Angle.m');
        keyboard
end


% Step 2: Find the number of cm on the screen

linear_size_cm = pixel_size / pixels_per_cm;

half_linear_size = linear_size_cm / 2;


% Step 3: Find the degrees subtended by the given linear size

half_angular_size_rad = atan( half_linear_size / monitor_distance );

half_angular_size_deg = (360/(2*pi)) * half_angular_size_rad;

angular_size = 2 * half_angular_size_deg;





