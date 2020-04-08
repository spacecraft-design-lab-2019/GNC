function B_NED = get_magnetic_field(lat, lon, alt, year, order)

%     /*
%     lat is geocentric latitude in degrees
%     lon is longitude in degrees
%     alt is altitude in km
%     year is the fractional year (include months/days essentially)
% 
%     outputs B vector in NED
%     */

    deg2rad = pi / 180;
    lat = 90 - lat;
    lat = lat*deg2rad;
    lon = lon*deg2rad;
%     // radius of earth
    a = 6371.2;
    r = a + alt;
%     // year since 2015 for secular variation
    dt = year - 2015.0;
%     // magnetic field components
    B_r = 0; 
    B_lat = 0;
    B_lon = 0;

%     // IGRF model B field calculation, n/m sets order (5)

% get coefficients
g = get_g_coefficients();
h = get_h_coefficients();
g_sv = get_g_sv_coefficients();
h_sv = get_h_sv_coefficients();
P = get_P_coefficients(cos(lat),order);
Pd = get_Pd_coefficients(P, cos(lat), order);

    for n = 1:order+1
        for m = 0:order+1
            coef = power((a/r), n + 2) * ((g(n+1, m+1) + dt * g_sv(n+1, m+1)) * cos(m * lon) +...
                    (h(n+1, m+1) + dt * h_sv(n+1, m+1)) * sin(m * lon));
%             // Radial component
            B_r =  B_r + (n + 1) * coef * P(n+1, m+1);
%             // Colatitudinal component
            B_lat = B_lat - coef * Pd(n+1, m+1);
%             // Address singularity at colatitude of 0
            if sin(lat) == 0
%                 // Longitudinal component
                B_lon = B_lon + -cos(lat) * power((a/r), n + 2) * (- (g(n+1, m+1) + dt * g_sv(n+1, m+1)) * sin(m * lon) +...
                        (h(n+1, m+1) + dt * h_sv(n+1, m+1)) * cos(m * lon)) * Pd(n+1, m+1);
            else
                B_lon = B_lon + (-1/sin(lat)) * power((a/r), n + 2) * m * (-(g(n+1, m+1) + dt * g_sv(n+1, m+1)) * sin(m * lon) +...
                        (h(n+1, m+1) + dt * h_sv(n+1, m+1)) * cos(m * lon)) * P(n+1, m+1);
            end
        end
    end


%     // NED (North, East, Down) coordinate frame
    B_vec(1) = -B_lat;
    B_vec(2) = B_lon;
    B_vec(3) = -B_r;

B_NED = B_vec;

end



function g = get_g_coefficients()

    g = [0.0,     0.0,    0.0,    0.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0;
         -29442.0, -1501.0,    0.0,    0.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0;
          -2445.1,  3012.9, 1676.7,    0.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0;
           1350.7, -2352.3, 1225.6,  582.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0;
            907.6,   813.7,  120.4, -334.9,   70.4,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0;
           -232.6,   360.1,  192.4, -140.9, -157.5,  4.1,   0.0,   0.0,  0.0,   0.0,  0.0;
             70.0,    67.7,   72.7, -129.9,  -28.9, 13.2, -70.9,   0.0,  0.0,   0.0,  0.0;
             81.6,   -76.1,   -6.8,   51.8,   15.0,  9.4,  -2.8,   6.8,  0.0,   0.0,  0.0;
             24.2,     8.8,  -16.9,   -3.2,  -20.6, 13.4,  11.7, -15.9, -2.0,   0.0,  0.0;
             5.4,      8.8,    3.1,   -3.3,    0.7, -13.3, -0.1,   8.7,  9.1, -10.5,  0.0;
            -1.9,     -6.3,    0.1,    0.5,   -0.5,   1.8, -0.7,   2.1,  2.4,  -1.8, -3.6];

end

function h = get_h_coefficients()
    h = [0.0,     0.0,    0.0,    0.0,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0;
         0.0,  4797.1,    0.0,    0.0,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0;
         0.0, -2845.6, -641.9,    0.0,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0;
         0.0,  -115.3,  244.9, -538.4,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0;
         0.0,   283.3, -188.7,  180.9, -329.5,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0;
         0.0,    47.3,  197.0, -119.3,   16.0, 100.2,   0.0,  0.0,  0.0,  0.0,  0.0;
         0.0,   -20.8,   33.2,   58.9,  -66.7,   7.3,  62.6,  0.0,  0.0,  0.0,  0.0;
         0.0,   -54.1,  -19.5,    5.7,   24.4,   3.4, -27.4, -2.2,  0.0,  0.0,  0.0;
         0.0,    10.1,  -18.3,   13.3,  -14.6,  16.2,   5.7, -9.1,  2.1,  0.0,  0.0;
         0.0,   -21.6,   10.8,   11.8,   -6.8,  -6.9,   7.8,  1.0, -4.0,  8.4,  0.0;
         0.0,     3.2,   -0.4,    4.6,    4.4,  -7.9,  -0.6, -4.2, -2.8, -1.2, -8.7];

end

function g_sv = get_g_sv_coefficients()

    g_sv = [ 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
            10.3,  18.1,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
            -8.7,  -3.3,  2.1,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             3.4,  -5.5, -0.7, -10.1,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
            -0.7,   0.2, -9.1,   4.1, -4.3,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
            -0.2,   0.5, -1.3,  -0.1,  1.4,  3.9,  0.0,  0.0, 0.0, 0.0, 0.0;
            -0.3,  -0.1, -0.7,   2.1, -1.2,  0.3,  1.6,  0.0, 0.0, 0.0, 0.0;
             0.3,  -0.2, -0.5,   1.3,  0.1, -0.6, -0.8,  0.2, 0.0, 0.0, 0.0;
             0.2,   0.0, -0.6,   0.5, -0.2,  0.4,  0.1, -0.4, 0.3, 0.0, 0.0;
             0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0];

end

function h_sv = get_h_sv_coefficients()

    h_sv = [ 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0, -26.6,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0, -27.4, -14.1,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0,   8.2,  -0.4,  1.8,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0,  -1.3,   5.3,  2.9, -5.2,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0,   0.6,   1.7, -1.2,  3.4,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0,   0.0,  -2.1, -0.7,  0.2,  0.9,  1.0,  0.0, 0.0, 0.0, 0.0;
             0.0,   0.8,   0.4, -0.2, -0.3, -0.6,  0.1, -0.2, 0.0, 0.0, 0.0;
             0.0,  -0.3,   0.3,  0.1,  0.5, -0.2, -0.3,  0.3, 0.0, 0.0, 0.0;
             0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;
             0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0];

end

function P = get_P_coefficients(x,order)

    P = zeros(order+2, order+2);
    P(1,1) = 1.0;
    P(2, 2) = sqrt(1 - power(x, 2));
    for n = 1:order+2
        for m = 0:order+2
            if (n ~= 1 || m ~= 1)
                prev2 = n - 2;
                if prev2 < 0
                    prev2 = 0;
                end
                if (m < n)
                    P(n+1, m+1) = (2 * n - 1) / sqrt(power(n, 2) - power(m, 2)) * x * P(n, m+1) - sqrt((power(n-1, 2) - power(m, 2)) / (power(n, 2) - power(m, 2))) * P(prev2+1, m+1);
                else
                    P(n+1, m+1) = sqrt(1.0 - 1.0/(2.0*m)) * sqrt(1 - power(x, 2)) * P(n , m );
                end
            end
        end
    end
end


function Pd = get_Pd_coefficients(P, x, order)

    Pd = zeros(order+1, order+1);
    Pd(2,2) = x;
    for n = 1:order+1
        for m = 0:order+1
            if (n~=1 || m ~= 1)
                if (m < n)
                    prev2 = n - 2;
                    if (prev2 < 0)
                        prev2 = 0;
                    end
                    Pd(n+1, m+1) = (2 * n - 1) / sqrt(power(n, 2) - power(m, 2)) * (x * Pd(n, m+1) - ...
                                                                            sqrt(1 - power(x, 2)) * P(n, m+1)) - sqrt((power(n-1, 2) - power(m, 2))/(power(n, 2) - power(m, 2))) * Pd(prev2+1, m+1);
                else
                    Pd(n+1, m+1) = sqrt(1.0 - 1.0 / (2.0 * m)) * (sqrt(1 - power(x, 2)) * Pd(n, m) + x * P(n , m));
                end
            end
        end
    end
end

