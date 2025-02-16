
clear,clc

dc = [38.907192, -77.036873];  %lat, lon
sf = [37.774929, -122.419418]; %lat, lon

re_mi = 3963 * sind(dc(1)) % mi

lat_diff_deg = dc(2) - sf(2)

dist_mi = deg2rad(lat_diff_deg) * re_mi

%%
clear, clc

% from maps.google.com
dc = [38.907192, -77.036873];  %lat, lon
sf = [37.774929, -122.419418]; %lat, lon

% radius of earth, mi
re_mi = 3963; 

% compute ECEF position vectors
dcECEF = computeECEF(dc(1), dc(2), re_mi);
sfECEF = computeECEF(sf(1), sf(2), re_mi);

% take dot product (not using the dot function since this is a hw assignment)
dc_dot_sc = dcECEF(1) * sfECEF(1) + dcECEF(2) * sfECEF(2) + dcECEF(3) * sfECEF(3);

% find angle, assuming magnitude of vectors are == r
theta = acos(dc_dot_sc / (re_mi * re_mi));

% get great circle distance (arclength)
dist = theta * re_mi


function ECEF = computeECEF(lat, lon, re)
% given lat, lon, and radius of earth, compute ECEF coordinate

    X = re * cosd(lat) * cosd(lon); 
    Y = re * cosd(lat) * sind(lon);
    Z = re * sind(lat);

    ECEF = [X, Y, Z];

end




