%% Generate timings and errors for a given problem size and given kernel
% parameters are defined in the paper 
% A Tutorial and Open Source Software for the Efficient Evaluation of Gravity and Magnetic Kernels (2019)
% Jarom Hogue, Rosemary Renaut and Saeed Vatankhah
% Trademarks: 
% Rosemary Renaut and Jarom Hogue (TM)
% Test_Efficiency Version 1: December 13, 2019. 
%%
nsx=25*scale;nsy=15*scale;nbz=2*scale;
padxl=padchoices(padding, 1);padxr=padchoices(padding, 2);padyl=padchoices(padding, 3);padyr=padchoices(padding, 4);
padxl=round(nsx*padxl/100);padxr=round(nsx*padxr/100);padyr=round(nsy*padyr/100);padyl=round(nsy*padyl/100);
gsx=80/scale;gsy=80/scale;gsz=100/scale;
Totalz=gsz*nbz;zOrigin=0;
Totalx=gsx*(nsx+padxr);xOrigin=-padxl*gsx;
Totaly=gsy*(nsy+padyr);yOrigin=-padyl*gsy;
nbx=(nsx+(padxr+padxl));                                    %  total number of blocks in x direction
nby=(nsy+(padyr+padyl));                                    %  total number of blocks in y direction
m   = nsx*nsy;                                              %  number of observation points
n   = nbx*nby*nbz;                                          %  number of model parameters
nr	= nbx*nby;                                              %  Column dimension block r
padx=padxl+padxr;                                           %  Total padding in x 
pady=padyl+padyr;                                           %  Total padding in y 
prob_params=[nsx,nsy,nbz,padxl,padxr,padyl,padyr,nbx,nby,m,n,nr,padx,pady];
z_blocks=gsz*(0:nbz);
%% Determine if matrix will be out of memory when calculating sensitivity matrix G
if scalemax==upperscale
    try zeros(m,n);
    catch ME
        if ME.identifier=='MATLAB:array:SizeLimitExceeded'
            scalemax=scale-1;
        end
    end
end
%% Find the optimal fft that is used in creating the matrices, forming Gx and G'x and save for future use 
rng('default');
fftw('dwisdom',[]);
switch prob_kind
    case 1
        fftw('planner','exhaustive');
        [That]=forward_gravity_bttb(gsx,gsy,z_blocks,prob_params);
        fftinfoG{scale}=fftw('dwisdom');
        fftw('dwisdom',[]);
        fftw('planner','exhaustive');
        Gx_bttb=matrix_mult_bttb(That,rand(n,1),1,prob_params);
        fftinfoGx{scale}=fftw('dwisdom');
        fftw('planner','exhaustive');
        Gx_bttb=matrix_mult_bttb(That,rand(m,1),2,prob_params);
        fftinfoGTx{scale}=fftw('dwisdom');
        fftw('dwisdom',[]);
        clear That
end
%% Run over totalruns for the given problem size and calculate the timings and differences in evaluating matrix-vector products
for nruns=1:totalruns
    mvec=rand(n,1);
    dvec=rand(m,1);
    switch prob_kind
        case 1
            fftw('dwisdom',fftinfoG{scale});
            tic;[That]=forward_gravity_bttb(gsx,gsy,z_blocks,prob_params);timing=toc;timeG_bttb(prob_kind,nruns,scale)=timing;
            fftw('dwisdom',fftinfoGx{scale});
            tic;true_data_bttb =matrix_mult_bttb(That,mvec,1,prob_params);timing=toc;timeGx_bttb(prob_kind,nruns,scale)=timing;
            fftw('dwisdom',fftinfoGTx{scale});
            tic;true_transpose_data_bttb =matrix_mult_bttb(That,dvec,2,prob_params);timing=toc;timeGTx_bttb(prob_kind,nruns,scale)=timing;
            clear That
            if scale<=scalemax
                tic;[G]=forward_gravity(gsx,gsy,z_blocks,prob_params);timing=toc;timeG(prob_kind,nruns,scale)=timing;
                tic;true_data =G*mvec;timing=toc;timeGx(prob_kind,nruns,scale)=timing;
                tic;true_transpose_data =G'*dvec;timing=toc;timeGTx(prob_kind,nruns,scale)=timing;
                clear G
                error_bttb(prob_kind,nruns,scale)=norm(true_data_bttb-true_data)/norm(true_data);
                error_transpose_bttb(prob_kind,nruns,scale)=norm(true_transpose_data_bttb-true_transpose_data)/norm(true_transpose_data);
            end
        case 2
            fftw('dwisdom',fftinfoG{scale});
            tic;[That]=forward_magnetic_bttb(gsx,gsy,z_blocks,prob_params,D,I,H);timing=toc;timeG_bttb(prob_kind,nruns,scale)=timing;
            fftw('dwisdom',fftinfoGx{scale});
            tic;true_data_bttb =matrix_mult_bttb(That,mvec,1,prob_params);timing=toc;timeGx_bttb(prob_kind,nruns,scale)=timing;
            fftw('dwisdom',fftinfoGTx{scale});
            tic;true_transpose_data_bttb =matrix_mult_bttb(That,dvec,2,prob_params);timing=toc;timeGTx_bttb(prob_kind,nruns,scale)=timing;
            clear That
            if scale<=scalemax
                tic;[G]=forward_magnetic(gsx,gsy,z_blocks,prob_params,D,I,H);timing=toc;timeG(prob_kind,nruns,scale)=timing;
                tic;true_data =G*mvec;timing=toc;timeGx(prob_kind,nruns,scale)=timing;
                tic;true_transpose_data =G'*dvec;timing=toc;timeGTx(prob_kind,nruns,scale)=timing;
                clear G
                error_bttb(prob_kind,nruns,scale)=norm(true_data_bttb-true_data)/norm(true_data);
                error_transpose_bttb(prob_kind,nruns,scale)=norm(true_transpose_data_bttb-true_transpose_data)/norm(true_transpose_data);
            end
    end
end