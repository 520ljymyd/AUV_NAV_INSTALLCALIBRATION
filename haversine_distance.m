function distance = haversine_distance(latlon1, latlon2)
% Function to calculate the distance between two points on Earth's surface
% Inputs:
%   - latlon1: [latitude, longitude] of point 1 (in degrees)
%   - latlon2: [latitude, longitude] of point 2 (in degrees)
% Output:
%   - distance: Distance in meters

    % Earth's radius in meters
    R = 6371000;

    % Convert latitude and longitude from degrees to radians
    latlon1 = deg2rad(latlon1);
    latlon2 = deg2rad(latlon2);

    % Haversine formula
    delta_lat = latlon2(1) - latlon1(1);
    delta_lon = latlon2(2) - latlon1(2);
    a = sin(delta_lat / 2)^2 + cos(latlon1(1)) * cos(latlon2(1)) * sin(delta_lon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distance = R * c;  % Distance in meters
end
