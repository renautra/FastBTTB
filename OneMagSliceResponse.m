function g = OneMagPrismResponse(z1,z2,X,Y,R,gc)
% This routine computes the response at one station for a given set of
% coordinates or  multiple stations for given coordinates. The dimensions
% are calculated inside the function. 
%
% The code is based on 
% Rao and Babu (1991) : A rapid method for three-dimensional modeling of
% magnetic anomalies in Geophysics
% For the geometry
%   ^ x (geographic North)
%   |
%   |
%   |--------------> y   (geographic East) z-axis points vertically downwards 
%
% Input Parameters are:
%%
% * $z_1, z_2$       depth coordinates of prism $z_1<z_2$. 
% * $X, Y$  arrays 
% * $R =  X(:)^2+Y^2$
%%
% Output g is the response for the given set of data in X,Y and R
%  
% Trademarks: 
%
% Rosemary Renaut 2019 (TM)
% Original code developer Saeed Vatankhah
%%
l=length(X)-1;k=length(Y)-1;
X=reshape(X,1,l+1);Y=reshape(Y,k+1,1);% to assure that the operations work
R1=sqrt(R+z1.^2);R2=sqrt(R+z2.^2);
F1=((R2(1:k,1:l)+ X(1:l))./(R1(1:k,1:l)+ X(1:l))).*((R1(1:k,2:l+1)+ X(2:l+1))./(R2(1:k,2:l+1)+ X(2:l+1))).*...
    ((R1(2:k+1,1:l)+ X(1:l))./(R2(2:k+1,1:l)+ X(1:l))).*((R2(2:k+1,2:l+1)+ X(2:l+1))./(R1(2:k+1,2:l+1)+ X(2:l+1)));     
F2=((R2(1:k,1:l)+ Y(1:k))./(R1(1:k,1:l)+ Y(1:k))).*((R1(1:k,2:l+1)+ Y(1:k))./(R2(1:k,2:l+1)+ Y(1:k))).*...
    ((R1(2:k+1,1:l)+ Y(2:k+1))./(R2(2:k+1,1:l)+ Y(2:k+1))).*((R2(2:k+1,2:l+1)+ Y(2:k+1))./(R1(2:k+1,2:l+1)+Y(2:k+1)));
F3=((R2(1:k,1:l)+ z2)./(R1(1:k,1:l)+ z1)).*((R1(1:k,2:l+1)+z1)./(R2(1:k,2:l+1)+z2)).*...
    ((R1(2:k+1,1:l)+ z1)./(R2(2:k+1,1:l)+z2)).*((R2(2:k+1,2:l+1)+ z2)./(R1(2:k+1,2:l+1)+ z1));
F4=atan2(X(2:l+1).*z2,R2(2:k+1,2:l+1).*Y(2:k+1))-atan2(X(1:l).*z2,R2(2:k+1,1:l).*Y(2:k+1))-atan2(X(2:l+1).*z2,R2(1:k,2:l+1).*Y(1:k))+...
    atan2(X(1:l).*z2,R2(1:k,1:l).*Y(1:k))-atan2(X(2:l+1).*z1,R1(2:k+1,2:l+1).*Y(2:k+1))+atan2(X(1:l).*z1,R1(2:k+1,1:l).*Y(2:k+1))+...
    atan2(X(2:l+1).*z1,R1(1:k,2:l+1).*Y(1:k))-atan2(X(1:l).*z1,R1(1:k,1:l).*Y(1:k));
F5=atan2(Y(2:k+1).*z2,R2(2:k+1,2:l+1).*X(2:l+1))-atan2(Y(2:k+1).*z2,R2(2:k+1,1:l).*X(1:l))-atan2(Y(1:k).*z2,R2(1:k,2:l+1).*X(2:l+1))+...
    atan2(Y(1:k).*z2,R2(1:k,1:l).*X(1:l))-atan2(Y(2:k+1).*z1,R1(2:k+1,2:l+1).*X(2:l+1))+atan2(Y(2:k+1).*z1,R1(2:k+1,1:l).*X(1:l))+...
    atan2(Y(1:k).*z1,R1(1:k,2:l+1).*X(2:l+1))-atan2(Y(1:k).*z1,R1(1:k,1:l).*X(1:l));           
g=(gc(1)*log(F1)+gc(2)*log(F2)+gc(3)*log(F3)+gc(4)*F4+gc(5)*F5);
g=g(:);