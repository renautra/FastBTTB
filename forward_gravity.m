function[G]=forward_gravity(gsx,gsy,z_blocks,prob_params)
% Generates the kernel matrix for the gravity problem.
% Based on Chen and Liu 2018
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
nX=nsx+max(padxl,padxr);
nY=nsy+max(padyl,padyr);
X=(-0.5:1:nX-0.5)*gsx; 
Y=(-0.5:1:nY-0.5)*gsy; 
X2=X.^2;Y2=Y.^2;
XY=X(:).*Y;
R=X2(:)+Y2;
G=zeros(m,n);
Gr=zeros(m,nr);
Grq=zeros(nsx,nbx,nY);
cb=1;
for hz=1:nbz
    g=OneGravSliceResponse(z_blocks(hz),z_blocks(hz+1),X,Y,XY,R);%gives a row of length(R)
    for j=1:nY
        r=g(((j-1)*nX+1):j*nX);
        Grq(:,:,j)=toeplitz(r(padxl+1:nsx+padxl), r([padxl+1:-1:2 1:nsx+padxr]));
    end
    cj=1; 
    for j=[padyl+1:-1:2 1:nsy+padyr]  Gr(1:nsx,cj:nbx+cj-1)= Grq(:,:,j);cj=cj+nbx;end
    i=nsx;
    for j=2:nsy
        Gr(i+1:i+nsx,:)=[Grq(:,:,j+padyl),Gr(i-nsx+1:i,1:(nby-1)*nbx)];
        i=i+nsx;
    end
    G(1:m,cb:cb+nr-1)=Gr;cb=cb+nr;
end