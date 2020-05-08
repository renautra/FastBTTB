function [That]=forward_gravity_bttb(gsx,gsy,z_blocks,prob_params)
% Generates the FFT matrix T for the gravity problem. Vogel 2002
% Based on Chen and Liu 2018
%%
% * gsx : grid size for x
% * gsy : grid size for y
% * z_blocks : coordinates for z grid
% * prob_param : parameters that define the problem
%
% Trademarks: 
% Jarom Hogue, October 2019(TM)
[nsx, nsy, nbz, padxl, padxr, padyl, padyr, nbx, nby]=matsplit(prob_params(1:9));
nX=nsx+max(padxl,padxr);
nY=nsy+max(padyl,padyr);
X=(-0.5:1:nX-0.5)*gsx; 
Y=(-0.5:1:nY-0.5)*gsy; 
X2=X.^2;Y2=Y.^2;
XY=X(:).*Y;
R=X2(:)+Y2;
Atemp=zeros(nsx+nbx-1,nsy+nby-1);
That=zeros(nsx+nbx-1,nsy+nby-1,nbz);
fftw('dwisdom',[]);             % Set up to find optimal fft
fftw('planner','exhaustive');   % Do exhaustive planning
for hz=1:nbz
    g=OneGravSliceResponse(z_blocks(hz),z_blocks(hz+1),X,Y,XY,R);
    ct=1;
    for j=[1+padyl:nsy+padyl nsy+padyr:-1:2 1:padyl]
        r=g(((j-1)*nX+1):j*nX);
        Atemp(:,ct) = [r(1+padxl:nsx+padxl);r(nsx+padxr:-1:2);r(1:padxl)];
        ct = ct+1;
    end
    if hz==1 
        That(:,:,hz) = fft2(Atemp);
        fftinfo=fftw('dwisdom');
        fftw('dwisdom',fftinfo);% each block uses the same size
    else
        That(:,:,hz) = fft2(Atemp);
    end
end
