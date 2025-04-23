% EN.525.645.82.SP25 Homework 8
% Brian Caskey 

%% 8.1.1 
% Reproduce the calculations and results of the sextant example in the 
% current Nautical Almanac. (See the PDF file attached above for guidance, 
% but use the current year data instead of that provided with the file.) 
% To check your work, the current Almanac and its example provides i the 
% solution.

% assumptions
initial_lat = deg2rad(-15); % W 
initial_lon = deg2rad(32);  % N 

% From the 2025 nautical almanac at:
% https://thenauticalalmanac.com/TNARegular/2025_Nautical_Almanac.pdf
%
% hour angles are for 2025 June 21. I interpolate linearly to get the
% lookup angles from the table

% ---- Regulus at 20h 29m 23s ---------------------------------------------
H0_regulus = deg2rad(37.4204);
t_hr = 29/60 + 23/3600; % how far through the hour the lookup is

GHA_aries_lookup(1) = deg2rad(210 + 160/60);    % 20hr value
GHA_aries_lookup(2) = deg2rad(225 + 18.5/60);   % 21hr value

% true GHA_aries value at lookup time
GHA_aries = GHA_aries_lookup(1) + (GHA_aries_lookup(2) - GHA_aries_lookup(1)) * t_hr;

% compute regulus GHA (wrap around 2pi) 
SHA_regulus = deg2rad(207 + 34.2/60);
GHA_regulus = GHA_aries + SHA_regulus - 2*pi;

% store Declination lookup
dec_regulus = deg2rad(11 + 50.6/60);

% ---- Antares at 20h 45m 47s ---------------------------------------------
H0_antares = deg2rad(20.3226);
t_hr = 45/60 + 47/3600; % how far through the hour the lookup is

GHA_aries_lookup(1) = deg2rad(210 + 160/60);    % 20hr value
GHA_aries_lookup(2) = deg2rad(225 + 18.5/60);   % 21hr value

% true GHA_aries value at lookup time
GHA_aries = GHA_aries_lookup(1) + (GHA_aries_lookup(2) - GHA_aries_lookup(1)) * t_hr;

% compute regulus GHA (wrap around 2pi) 
SHA_antares = deg2rad(112 + 15.0/60);
GHA_antares = GHA_aries + SHA_antares;

% store Declination lookup
dec_antares = deg2rad(-26 -29/60);

% ---- Kochab  at 21h 10m 34s ---------------------------------------------
H0_kochab = 47.2050;
t_hr = 10/60 + 34/3600; % how far through the hour the lookup is

GHA_aries_lookup(1) = deg2rad(225 + 18.5/60);    % 21hr value
GHA_aries_lookup(2) = deg2rad(240 + 21.0/60);   % 22hr value

% true GHA_aries value at lookup time
GHA_aries = GHA_aries_lookup(1) + (GHA_aries_lookup(2) - GHA_aries_lookup(1)) * t_hr;

% compute regulus GHA 
SHA_kochab = deg2rad(137 + 19.0/60);
GHA_kochab = GHA_aries + SHA_kochab;

% store Declination lookup
dec_kochab = deg2rad(74 + 03.2/60);

% =========================================================================

% store inputs in to array
GHA = [GHA_regulus, GHA_antares, GHA_kochab];
dec = [dec_regulus, dec_antares, dec_kochab];
H0  = [H0_regulus, H0_antares, H0_kochab];

lon = initial_lon;
lat = initial_lat;

d = 100; % dummy start value

for idx = 1:10

    % compute local hour angle
    LHA = GHA - lon;
    
    % Compute calculated altitude (Hc), azimuth (Z), intercept (p) for each sighting
    for i = 1:3
        Hc(i) = asin(sin(dec(i) * sin(lat) + cos(dec(i) * cos(lat) * cos(LHA(i)))));
        X = (sin(dec(i)) - sin(lat) * sin(Hc(i))) / (cos(lat) * cos(Hc(i))); 
        if X > 1
            X = 1;
        elseif X < -1
            X = -1;
        end
        Z(i)  = acos(X);
        p(i)  = H0(i) - Hc(i);
    end

    % apply vector summation
    A = sum(cos(Z).^2);
    B = sum(cos(Z).^2 .* sin(Z));
    C = sum(sin(Z).^2);
    D = sum(p .* cos(Z));
    E = sum(p .* sin(Z));

    % improved estiamte
    G = A*C - B^2;

    % I'm guessing these are in degrees? 
    dLat_deg = (C*D - B*E)/G;
    dLon_deg = (A*E - B*D)/(G*cos(lat));

    dLat = deg2rad(dLat_deg);
    dLon = deg2rad(dLon_deg);

    lat_new = lat + dLat;
    lon_new = lon + deg2rad(dLon);

    % check distance
    d = 60 * sqrt(dLon * cos(lat)^2 + (dLat)^2);

    fprintf("distance = %f\n", d);

    % update lat & lon
    lat = lat_new;
    lon = lon_new;

    if d < 20
        break
    end

end

%% 8.1.3
% Calculate the Hc and azimuth to each star for this actual location and 
% time. Compare, in tabular form, to your estimates from above. These 
% calculated values of Hc are to be used as the values of Ho in your 
% computations

clear,clc

% My position
myLat = deg2rad(39.284017);
myLon = deg2rad(-76.610119);

% get Aries GHA data from almanac
GHA_aries(1) = deg2rad(209 + 19.3/60); % GHA at 00h 2025.04.21
GHA_aries(2) = deg2rad(224 + 21.7/60); % GHA at 01h 2025.04.21

% ---- Sirius -------------------------------------------------------------
SHA_sirius = deg2rad(258 + 26.2/60); % from almanac 2025.04.21
dec_sirius = deg2rad(-16 - 45.2/60); % from almanac 2025.04.21

tHr = 21/60 + 31/3600; % how far through the hour the sighting is

% compute GHA sirius as SHA + interpolated GHA_aries value
GHA_sirius = SHA_sirius + GHA_aries(1) + (GHA_aries(2) - GHA_aries(1))*tHr;

% compute local hour angle (rad)
LHA_sirius = GHA_sirius + myLon - 2*pi;

% compute S, C, and altitude Hc from
S = sin(dec_sirius);
C = cos(dec_sirius) * cos(LHA_sirius);
Hc_sirius = asin(S * sin(myLat) + C * cos(myLat));

% Compute azimuth
X = (S * cos(myLat) - C * sin(myLat)) / cos(Hc_sirius);
A = acos(X);

if LHA_sirius > pi
    Z_sirius = A;
else
    Z_sirius = 2*pi-A;
end

fprintf('Sirius altitude (deg) = %f\n', rad2deg(Hc_sirius));
fprintf('Sirius azimuth (deg)  = %f\n', rad2deg(Z_sirius));
fprintf('-----------------------------------------\n');

% ---- Arcturus -----------------------------------------------------------
SHA_arcturus = deg2rad(145 + 47.4/60); % from almanac 2025.04.21
dec_arcturus = deg2rad(019 + 02.9/60); % from almanac 2025.04.21

tHr = 22/60 + 09/3600; % how far through the hour the sighting is

% compute GHA sirius as SHA + interpolated GHA_aries value
GHA_arcturus = SHA_arcturus + GHA_aries(1) + (GHA_aries(2) - GHA_aries(1))*tHr;

% compute local hour angle (rad)
LHA_arcturus = GHA_arcturus + myLon - 2*pi;

% compute S, C, and altitude Hc from
S = sin(dec_arcturus);
C = cos(dec_arcturus) * cos(LHA_arcturus);
Hc_arcturus = asin(S * sin(myLat) + C * cos(myLat));

% Compute azimuth
X = (S * cos(myLat) - C * sin(myLat)) / cos(Hc_arcturus);
A = acos(X);

if LHA_sirius > pi
    Z_arcturus = A;
else
    Z_arcturus = 2*pi-A;
end

fprintf('Arcturus altitude (deg) = %f\n', rad2deg(Hc_arcturus));
fprintf('Arcturus azimuth (deg)  = %f\n', rad2deg(Z_arcturus));
fprintf('-----------------------------------------\n');


% ---- Capella ------------------------------------------------------------
SHA_capella = deg2rad(280 + 21.9/60); % from almanac 2025.04.21
dec_capella = deg2rad(046 + 01.5/60); % from almanac 2025.04.21

tHr = 23/60 + 55/3600; % how far through the hour the sighting is

% compute GHA sirius as SHA + interpolated GHA_aries value
GHA_capella = SHA_capella + GHA_aries(1) + (GHA_aries(2) - GHA_aries(1))*tHr;

% compute local hour angle (rad)
LHA_capella = GHA_capella + myLon - 2*pi;

% compute S, C, and altitude Hc from
S = sin(dec_capella);
C = cos(dec_capella) * cos(LHA_capella);
Hc_capella = asin(S * sin(myLat) + C * cos(myLat));

% Compute azimuth
X = (S * cos(myLat) - C * sin(myLat)) / cos(Hc_capella);
A = acos(X);

if LHA_sirius > pi
    Z_capella = A;
else
    Z_capella = 2*pi-A;
end

fprintf('Capella altitude (deg) = %f\n', rad2deg(Hc_capella));
fprintf('Capella azimuth (deg)  = %f\n', rad2deg(Z_capella));
fprintf('-----------------------------------------\n');


%% 8.1.4
% 4 Calculate new values of Hc for each star under the assumption that you
% are at 0 degrees longitude and 0 degrees latitude for the same instant in
% time. These incorrect values are then used, with your Ho values, to drive
% the algorithm in the Nautical Almanac (from the PDF file referenced 
% above). Using several iterations of this algorithm, show that the 
% algorithm eventually converges to your correct position.  

% NOTE THE SECION BEFORE THIS MUST BE RUN FOR THIS CODE TO WORK

% store inputs from section above into array
GHA = [GHA_sirius, GHA_arcturus, GHA_capella];
dec = [dec_sirius, dec_arcturus, dec_capella];
H0  = [Hc_sirius, Hc_arcturus, Hc_capella];

% initial guesses
lon = 0;
lat = 0;

d = 100; % dummy start value

for idx = 1:100

    % compute local hour angle
    LHA = GHA + lon;
    
    % Compute calculated altitude (Hc), azimuth (Z), intercept (p) for each sighting
    for i = 1:3
        S = sin(dec(i));
        C = cos(dec(i)) * cos(LHA(i));
        Hc(i) = asin(S * sin(myLat) + C * cos(myLat));
        X = (S * cos(myLat) - C * sin(myLat)) / cos(Hc(i)); 
        if X > 1
            X = 1;
        elseif X < -1
            X = -1;
        end
        z(i)  = acos(X);
        p(i)  = rad2deg(H0(i)) - rad2deg(Hc(i)); % nm
    end

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

%% 8.3.1
% In navigation problems, it is common to divide by the sine or cosine of 
% an angle. For a 1% measurement error in the measurement of θ, plot the 
% corresponding estimation error in the function f(θ) = 1/cos(θ) 
% for 0 < θ < 90 degrees.

clear, clc

% plot range
theta = 0:0.001:pi; % rads

% 1% measurement error 
theta_error = 0.01 .* theta;

% compute error vals
f_nominal = 1 ./ cos(theta);
f_plus = 1 ./ cos(theta + theta_error);
f_minus = 1 ./ cos(theta - theta_error);

error_plus  = f_plus - f_nominal;
error_minus = f_minus - f_nominal;


% plot
figure; hold on; grid on;
lw = 1.5;
plot(rad2deg(theta), abs(error_plus), '-k', 'LineWidth', lw); lgnd(1) = "+1% Error";
plot(rad2deg(theta), abs(error_minus), '--k', 'LineWidth', lw); lgnd(end+1) = "-1% Error";

yscale log

xlabel('\theta (degrees)');
ylabel('Estimation Error in f(\theta)');
title('Estimation Error in f(\theta) = 1/cos(\theta) due to 1% Error in \theta');
grid on;
legend(lgnd)

