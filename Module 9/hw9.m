% 
clear, clc

% initial guesses
lat = deg2rad(40);
lon = deg2rad(-76);

% sighting values (computed from previous problem)
%             27W       115W
az = deg2rad([118.03, 360-128.80]);
H0 = deg2rad([21.93,  30.07]);

d = 100; % dummy start value

for idx = 1:100

    [~, Hc(1), ~, Hc(2)] = computeAzEl(lat, lon);
    p  = rad2deg(H0) - rad2deg(Hc); % nm
    z = az;

    % apply vector summation
    A = sum(cos(z).^2);
    B = sum(cos(z) .* sin(z));
    C = sum(sin(z).^2);
    D = sum(p .* cos(z));
    E = sum(p .* sin(z));

    % improved estiamte
    G = A*C - B^2;

    % I'm guessing these are in degrees since p is in nm?  
    dLat_deg = (C*D - B*E)/G;
    dLon_deg = (A*E - B*D)/(G*cos(lat));

    dLat = deg2rad(dLat_deg);
    dLon = deg2rad(dLon_deg);

    lat_new = lat + dLat;
    lon_new = lon + dLon;

    % check distance
    d = 60 * sqrt((dLat_deg)^2 * cos(lat)^2 + (dLat_deg)^2);

    % update lat & lon
    lat = lat_new;
    lon = lon_new;

    % bound lat & lon within 2pi
    if lat > 2*pi
        lat = lat-2*pi;
    end

    if lon > 2*pi
        lon = lon-2*pi;
    end

    fprintf('Step #%f\n', idx);
    fprintf("distance = %f\n", d);
    fprintf('Lat: %f\n', rad2deg(lat));
    fprintf('Lon: %f\n', rad2deg(lon));
    fprintf('------------------\n');

    if d < 20
        break
    end
end

fprintf('------------------\n');
fprintf('Final Lat: %f\n', rad2deg(lat));
fprintf('Final Lon: %f\n', rad2deg(lon));


function [az1, el1, az2, el2] = computeAzEl(myLat, myLon)
    
    %---- Fist satellite at 27W -----------------------------------------------
    satLon = deg2rad(-27); % degrees west
    
    lonDiff = satLon - myLon;
    
    Re = 6378;   % Earth radius, km
    Rs = 42164;  % Geostationary satellite orbit radius, km
    
    eta_s = asin( cos(myLat) * cos(lonDiff) );
    
    el1 = atan( (sin(eta_s) - (Re/Rs)) / cos(eta_s) );
    
    az1 = atan2( sin(lonDiff), -sin(myLat) * cos(lonDiff) );
    
    %---- Second satellite at 115W --------------------------------------------
    satLon = deg2rad(-115); % degrees west

    lonDiff = satLon - myLon;

    eta_s = asin( cos(myLat) * cos(lonDiff) );

    el2 = atan( (sin(eta_s) - (Re/Rs)) / cos(eta_s) );

    az2 = atan2( sin(lonDiff), -sin(myLat) * cos(lonDiff) );

end
