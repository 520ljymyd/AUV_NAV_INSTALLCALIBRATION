function [x, y, z] = latlon2xyz(latlon, R)
    % Convert latitude and longitude to Cartesian coordinates
    % Inputs:
    % - latlon: [latitude, longitude] in degrees
    % - R: radius of the Earth (default 6371 km for Earth)
    % Outputs:
    % - x, y, z: Cartesian coordinates

    lat = deg2rad(latlon(1));
    lon = deg2rad(latlon(2));

    x = R * cos(lat) * cos(lon);
    y = R * cos(lat) * sin(lon);
    z = R * sin(lat);
end