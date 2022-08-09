function [ dx ] = clino_dyn( t, x, const )

%  State differential
dx      = zeros(4,1);
dx(1:2) = x(3:4);
dx(3)   = -const.f/const.m*x(3) + const.mstar/const.m*const.omega^2*x(1) - const.mstar/const.m*const.g*sin(const.omega*t);
dx(4)   = -const.f/const.m*x(4) + const.mstar/const.m*const.omega^2*x(2) - const.mstar/const.m*const.g*cos(const.omega*t);

end