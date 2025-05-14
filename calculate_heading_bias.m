function [delta_OAB, angle_deg] = calculate_heading_bias(positionO, positionDR, positionTrue)
% Calculate heading bias based on DR and true trajectory positions
%
% Inputs:
%   - positionO: Reference point O [latitude, longitude] in degrees
%   - positionDR: DR trajectory end point [latitude, longitude] in degrees
%   - positionTrue: True trajectory end point [latitude, longitude] in degrees
% Outputs:
%   - delta_OAB: Scale factor difference between OA and OB
%   - angle_deg: Heading bias angle in degrees

    % Calculate distances using Haversine formula
    distanceOA = haversine_distance(positionO, positionDR);
    distanceOB = haversine_distance(positionO, positionTrue);
    distanceAB = haversine_distance(positionDR, positionTrue);

    % Scale factor difference
    delta_OAB = distanceOA / distanceOB - 1;

    % Calculate the heading bias angle using cosine rule
    cos_angle = (distanceOA^2 + distanceOB^2 - distanceAB^2) / (2 * distanceOA * distanceOB);
    cos_angle = max(min(cos_angle, 1), -1); % Ensure within [-1, 1] range for acos
    angle_rad = acos(cos_angle);  % Angle in radians
    angle_deg = rad2deg(angle_rad);  % Convert to degrees
end

function distance = haversine_distance(latlon1, latlon2)
% Calculate the Haversine distance between two points on Earth's surface
%
% Inputs:
%   - latlon1: [latitude, longitude] of point 1 in degrees
%   - latlon2: [latitude, longitude] of point 2 in degrees
% Output:
%   - distance: Distance in meters

    % Earth's radius in meters
    R = 6371000;

    % Convert degrees to radians
    latlon1 = deg2rad(latlon1);
    latlon2 = deg2rad(latlon2);

    % Haversine formula
    delta_lat = latlon2(1) - latlon1(1);
    delta_lon = latlon2(2) - latlon1(2);
    a = sin(delta_lat / 2)^2 + cos(latlon1(1)) * cos(latlon2(1)) * sin(delta_lon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distance = R * c;  % Distance in meters
end
