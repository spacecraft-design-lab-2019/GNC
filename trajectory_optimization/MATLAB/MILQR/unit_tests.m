clear 
clc
Earth = InitializeEarth();

% test custom IGRF function
    addpath('igrf')

    time = '08-Apr-2015 13:37:14';

    lat = 0;
    lon = 0;
    alt = 500;
    coord = 'geocentric';

    [Bx, By, Bz] = igrf(time, lat, lon, alt+Earth.r, coord);
    B_ned = get_magnetic_field(lat,lon,alt,2015,5);

    assert(B_ned(1)/Bx < 1.01 || B_ned(1)/Bx > .99 );
    assert(B_ned(2)/By < 1.01 || B_ned(2)/By > .99 );
    assert(B_ned(3)/Bz < 1.01 || B_ned(3)/Bz > .99 );
