function g = OneGravSliceResponse(z1,z2,X,Y,XY,R)
% This routine computes the response at one station for a given slice
% The code is based on 
% Fast and accurate forward modelling of gravity field using prismatic grids
% Longwei Chen and Lanbo Liu, (GJI - 2019) 
% For the geometry
%   ^ x (geographic East)
%   |
%   |
%   |--------------> y   (geographic North) z-axis points vertically downwards 
%
% Input Parameters are:
%%
% * $z_1, z_2$       depth coordinates of prism $z_1<z_2$. 
% * $X, Y$ are as stated in Chen and Liu. Created for all slices. 
% * $XY$ is an elementwise product $X(:)*Y$
% * $R =  X(:)^2+Y^2$
%%
% Output g is the response for a given slice
%  
% Trademarks: 
%
% Rosemary Renaut 2019 (TM)
% Original code developer Saeed Vatankhah
gamma=6.673844e-3; %Gravitational constant
[nx,ny]=size(R);
R1=sqrt(R+z1^2);
R2=sqrt(R+z2^2);
CMX=(log((X(:)+R1)./(X(:)+R2))).*Y;
CMY=(log((Y+R1)./(Y+R2))).*X(:);
CM5Z=atan2(XY,R1.*z1).*z1;
CM6Z=atan2(XY,R2.*z2).*z2;
CM56=CM5Z-CM6Z;
CM=(CM56-CMY-CMX)*gamma;
g=-(CM(1:nx-1,1:ny-1)-CM(1:nx-1,2:ny)-CM(2:nx,1:ny-1)+CM(2:nx,2:ny));
g=g(:);
