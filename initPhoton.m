function [x,y,z,ux,uy,uz] = initPhoton(type,ang)
% this is a function that initialize photons
if type == 1
    x = 0;y = 0;z = 0;
    ux = sind(ang);uy = 0;uz = cosd(ang);
end
end