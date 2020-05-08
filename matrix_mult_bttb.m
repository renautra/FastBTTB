function b=matrix_mult_bttb(That,x,t,prob_params)
%%
% * That          : Array - third dimension for each block in space 
% * x             : right hand side data
% * t             : 1: matrix multiplication, 2: transpose multiplication 3: x^TA
% * prob_params   : Parameters of the problem
% * b             : Output Gx,  G'Tx, or x'G
%
% Trademarks:
% Jarom Hogue, October 2019(TM).
%%
%global infoforward infoinverse
[nsx, nsy, nbz, nbx, nby ,m, n, nr, px, py]=matsplit(prob_params([1:3 8:14]));
W=zeros(nbx+nsx-1,nby+nsy-1);
if t>1
    W(1:nsx,1:nsy)=reshape(x,nsx,nsy);
    W=fft2(W);
    b=zeros(n,1);
else
    b=zeros(m,1);
end
for j = 1:nbz
    switch t
        case 1 %Gx
            W(1:nbx,1:nby)=reshape(x((j-1)*nr+1:j*nr),nbx,nby);
            W=real(ifft2(That(:,:,j).*fft2(W)));
            b=b+reshape(W(1:nsx,1:nsy),m,1);
            W=zeros(nbx+nsx-1,nby+nsy-1);
        case 2 %G'x
            Z=real(ifft2(conj(That(:,:,j)).*W));
            b((j-1)*nr+1:j*nr)=reshape(Z(1:nbx,1:nby),nr,1);
        case 3 %x'G
            Z=real(ifft2(conj(That(:,:,j)).*W));
            b((j-1)*nr+1:j*nr)=reshape(Z(1:nbx,1:nby),1,nr);
    end
end
