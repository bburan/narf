% *****************************************************************************
% Description: 
% Converts a GPS latitude or longitude degrees minutes string to a decimal
% degrees latitude or longitude
% *****************************************************************************
%
% Returns -1 if called with input argument in wrong format
% Returns -2 if not called with input arguments
% 
% Inputs:       degmin - latitude or longitude string to convert
%               EWNS - East, West, North, South indicator from NMEA string,
%               or use a NULL field ('') if not used
% Returns:      dec - decimal degrees version of input
% 
% *****************************************************************************
% Revision History:
%  1 - Initial release
%  2 - Use of E,W,N & S to indicate a positive or negative latitude or
%  longitude added
%  3 - Number of decimal places of output increased from default 4 to a max of 10 
%      (due to num2str function in line 53).
% *****************************************************************************
% Version: 3
% Release Date: 22/05/2003
% Author: Steve Dodds
% E-Mail: stevedodds@dsl.pipex.com
% *****************************************************************************
% Function convert latitude and longitude from degrees minutes to decimal degrees number
function dec = DegMin_2_Dec(degmin,EWNS)
% Latitude string format: ddmm.mmmm (dd = degrees)
% Longitude string format: dddmm.mmmm (ddd = degrees)

if nargin ~= 2
    dec = '-2'
else
    % Determine  if data is latitude or longitude
    switch length(strtok(degmin,'.'))
        case 4
            % latitude data
            deg = str2num(degmin(1:2)); % extract degrees portion of latitude string and convert to number
            min_start = 3;              % position in string for start of minutes
        case 5
            % longitude data
            deg = str2num(degmin(1:3)); % extract degrees portion of longitude string and convert to number
            min_start = 4;              % position in string for start of minutes
        otherwise
            % data not in correct format
            dec = '-1';
            return;
    end

    minutes = (str2num(degmin(min_start:length(degmin))))/60; % convert minutes to decimal degrees

    dec = num2str(deg + minutes,'%11.10g'); % degrees as decimal number

    % add negative sign to decimal degrees if south of equator or west
    switch EWNS
        case 'S'
            % south of equator is negative, add -ve sign
            dec = strcat('-',dec);
        case 'W'
            % west of Greenwich meridian is negative, add -ve sign
            dec = strcat('-',dec);
        otherwise
            % do nothing
    end
end
