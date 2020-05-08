function[G]=forward_magnetic(gsx,gsy,z_blocks,prob_params,D,I,H)
% Generates the kernel matrix for the magnetic problem.  
% Based on Rao and Babu (1991): Model parameters are inputs declination, inclination and scaled intensity of magnetization 
%%
% $D$, $I$, $H$
%%
% Note that the constants  are :
%%
% 
% $$\mu_0=4\pi 1e-7,\quad  C_m=\frac{\mu_0}{4*pi},\quad    H_1=1e-9\frac{F}{\mu_0},\quad    H=C_m 1e9 H_1=\frac{F}{4\pi}. $$
% 
% 
%
% Here $F$ is the intensity of the geomagnetic field in nT. $C_m$ in SI
% units, $\mu_0$ the magnetic permeability of free space. $H_1$ is magnetic
% field intensity (A/m). $H$ is the scaled magnetic field intensity. 
% Constants [G1 G2 G3 G4 G5] of Rao and Babu are multiplied by H and calculated one time only. Stored as [g1 g2 g3 g4 g5]
%%
% * gsx : grid size for x
% * gsy : grid size for y
% * z_blocks : coordinates for z grid
% * prob_param : parameters that define the problem
%
% Trademarks: 
% Rosemary Renaut 2019 (TM)
% Original code developer Saeed Vatankhah
[nsx, nsy, nbz, padxl, padxr, padyl, padyr, nbx, nby, m, n, nr]=matsplit(prob_params(1:12));
D=deg2rad(D); I=deg2rad(I);
I0=I; D0=D;
L=cos(I0)*cos(D0);
M=cos(I0)*sin(D0);
N=sin(I0);
g1=(2*M*N)*H;                  
g2=(2*L*N)*H;
g3=(2*M*L)*H;
g4=(N+M)*(N-M)*H; 
g5=(N+L)*(N-L)*H; 
gc=[g1,g2,g3,g4,g5];
nX=nsx+max(padxl,padxr);
nY=nsy+max(padyl,padyr);
X=(-(nX-0.5):1:nX-0.5)*gsx; 
Y=(-(nY-0.5):1:nY-0.5)*gsy;  
Y2=Y.^2;X2=X.^2;
R=X2(:)+Y2;
G=zeros(m,n);
Gr=zeros(m,nr);
grow = cell(nsy+padyl);
grow{nsy+padyl}=[];
gcol = cell(nsy+padyr);
gcol{nsy+padyl}=[];
cb=1;
warning('off')% note that due to Floating Point arithmetic the Toeplitz entries row and col position 1 are not exactly same
for hz=1:nbz
    ly=nY;
    grow{1} = OneMagSliceResponse(z_blocks(hz),z_blocks(hz+1),Y(ly:ly+nY),X(nX:2*nX)',R(nX:2*nX,ly:ly+nY),gc);
    for j=2:nsy+padyl 
        ly=ly-1; 
        grow{j} = OneMagSliceResponse(z_blocks(hz),z_blocks(hz+1),Y(ly:ly+1),X(nX:nX+nsx+padxr)',R(nX:nX+nsx+padxr, ly:ly+1),gc);
    end
    uy=nY+1;
    gcol{1} = OneMagSliceResponse(z_blocks(hz),z_blocks(hz+1),Y(uy:-1:uy-nY),X(nX+1:-1:1)',R(nX+1:-1:1,uy:-1:uy-nY),gc);
    for q=2:nsy+padyr 
        uy=uy+1; 
        gcol{q} = OneMagSliceResponse(z_blocks(hz),z_blocks(hz+1),Y(uy:-1:uy-1),X(nX+1:-1:1)',R(nX+1:-1:1, uy:-1:uy-1),gc);
    end
    cj=1; 
    for j=padyl+1:-1:2
        Gr(1:nsx,cj:nbx+cj-1)= toeplitz(gcol{1}((j-1)*nX+1+padxl:((j-1)*nX +nsx+padxl)),[gcol{1}(((j-1)*nX+1+padxl):-1:(j-1)*nX+2);grow{j}(1:nsx+padxr)]);
        cj=cj+nbx;
    end
    for q=1:nsy+padyr
         Gr(1:nsx,cj:nbx+cj-1)=toeplitz(gcol{q}(1+padxl:nsx+padxl),[gcol{q}(1+padxl:-1:2);  grow{1}((q-1)*nX+1:(q-1)*nX+nsx+padxr)]);
         cj=cj+nbx;
    end
    i=nsx;
    for j=padyl+2:nsy+padyl  
        Gr(i+1:i+nsx,:)=[toeplitz(gcol{1}((j-1)*nX+1+padxl:((j-1)*nX +nsx+padxl)),[gcol{1}((j-1)*nX+1+padxl:-1:(j-1)*nX+2);grow{j}(1:nsx+padxr)]),Gr(i-nsx+1:i,1:end-nbx)];
        i=i+nsx;
    end
    G(1:m,cb:cb+nr-1)=Gr;cb=cb+nr;
end
warning('on')