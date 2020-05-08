%% This script can be used to generate the results presented in the paper 
% A Tutorial and Open Source Software for the Efficient Evaluation of Gravity and Magnetic Kernels (2019)
% Jarom Hogue, Rosemary Renaut and Saeed Vatankhah
% Results presented in this paper were obtained using Matlab Release 2019b implemented on a computer wit Intel(R) Xeon (R) Gold $6138 $ processor  (2.00 GHz) and   256 GB RAM.
% The timing obtained on a different computer may vary, the differences (errors) may also be computer dependent. 
% Further there is a test builtin to prevent running with matrices that are too large. 
% For the noted computer we could run the full matrices up % to largest scale of 7. 
% On an iMac with 4.2Ghz Quad Core Intel Core i7 and 32GB RAM the largest scale is 6.
% Additional codes required to run this script are provided in the package.
% Trademarks: 
% Rosemary Renaut and Jarom Hogue (TM)
% TestingScript Version 1: December 13, 2019.  
%% Parameters that define the tests
% Users can change the parameters
% lowerscale, upperscale:     the smallest and largest values of scale used in the problem size [nsx nsy nbz]=[25 15 2]*scale
% scalemax:                   set to upperscale, but can be picked to restrict largest size of G directly generated
% padding:                    the percentage of padding according to padchoices
% totalruns:                  number of runs for the generation of the mean results   
%%
close all, clear all,
D=2;I=50;F=47000;H=(F)/(4*pi); %  H is scaled magnetic intensity (A/m)
%% Initialize parameters for the tests
upperscale=2;lowerscale=1;scalemax=upperscale;
padding=1;                     %for picking a padding from choices
padchoices=[0 0 0 0 ; 5 5 5 5; 2 2 0 0 ;0 0 2 2; 2 2 2 2;1 1 2 2; 2 2 1 1;1 2 1 1; 1 1 1 2;2 1 1 1; 2 1 2 1; 1 2 1 2];
totalruns=100;                 % Number of runs for finding mean values
timeG=zeros(2,totalruns,upperscale-lowerscale+1);timeGx=timeG;timeGTx=timeG;
timeG_bttb=timeG;timeGx_bttb=timeG;timeGTx_bttb=timeG;
error_bttb=timeG;error_transpose_bttb=timeG;
%% Run the tests for all problem sizes and both kernels
for prob_kind=1:2 %1 for gravity kernel and 2 for magnetic kernel
    for scale=lowerscale:upperscale
        Test_Efficiency;
    end
end
%% Generate results as tables for the screen - examine results online and generate data for plotting
S.type = '()';
rowct=1;
for prob_kind=1:2
    for scale=lowerscale:upperscale
        S.subs={prob_kind, ':',scale};        
        padxl=padchoices(padding, 1);padxr=padchoices(padding, 2);padyl=padchoices(padding, 3);padyr=padchoices(padding, 4);
        nsx=25*scale;nsy=15*scale;nbz=2*scale;
        padxl=round(nsx*padxl/100);padxr=round(nsx*padxr/100);padyr=round(nsy*padyr/100);padyl=round(nsy*padyl/100);
        col1table(rowct,:)=[prob_kind,(nsx),(nsy),2*scale,padxl,padxr,padyl,padyr,nsx*nsy,(nsx+(padxr+padxl))*(nsy+(padyr+padyl))*nbz];
        rowmeans(rowct,:)=[mean(subsref(timeG,S)),mean(subsref(timeG_bttb,S)),mean(subsref(timeGx,S)),mean(subsref(timeGTx,S)),mean(subsref(timeGx_bttb,S)),mean(subsref(timeGTx_bttb,S))];
        rowstds(rowct,:)=[std(subsref(timeG,S)),std(subsref(timeG_bttb,S)),std(subsref(timeGx,S)),std(subsref(timeGTx,S)),std(subsref(timeGx_bttb,S)),std(subsref(timeGTx_bttb,S))];
        errormean(rowct,:)=[mean(subsref(error_bttb,S)), mean(subsref(error_transpose_bttb,S))];
        errorstd(rowct,:)=[std(subsref(error_bttb,S)), std(subsref(error_transpose_bttb,S))];
        rowct=rowct+1;
    end
end
Table_Variable_Names={'Problem','sx','sy','nz','padxl','padxr','padyl','padyr','m','n','Cost_G','Cost_Gfft','Cost_Gx','Cost_GTx','Cost_Gfftx','Cost_GTfftx','Mean_Error_Forward','Mean_Error_Transpose','Std_Error_Forward','Std_Error_Transpose'};
T_mean_times=array2table([col1table, rowmeans],'VariableNames',{Table_Variable_Names{1:16}})
T_std_times=array2table([col1table, rowstds],'VariableNames',{Table_Variable_Names{1:16}})
T_errors=array2table([col1table, errormean, errorstd],'VariableNames',{Table_Variable_Names{1:10},Table_Variable_Names{17:20}})
%% Save the results for future upload to plot the results
save(['Results',int2str(padding)],'error_bttb','error_transpose_bttb','timeG_bttb','timeGx_bttb','timeGTx_bttb','timeG','timeGx','timeGTx','padding','padchoices',...
   'totalruns','T_mean_times','T_std_times','T_errors','col1table','rowmeans','rowstds','errormean','errorstd','upperscale','lowerscale','Table_Variable_Names','scalemax')