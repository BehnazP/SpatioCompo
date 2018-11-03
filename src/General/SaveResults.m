function SaveResults(times,models,pathname,burnin,skip,flag_save_csv)
% SAVERESULTS save the results into .mat file and .csv format
%
% SaveResults(times,models,pathname,burnin,skip,flag_save_csv)
% input:
% times             number of times to be saved
% models            different models to be saved
% pathname          the path where the results should be saved
% burnin            burn in samples of MCMC
% skip              skip samples in MCMC 
% flag_save_csv     indication of saving the files in .csv format;default 0
%
% SAVERESULTS.m 2018-07-15 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

if nargin < 3
    burnin = 1e4;
    skip = 10;
    flag_save_csv = 0;
end
    
if isempty(burnin)
    burnin = 1e4;
    skip = 10;
end

if isempty(models) && isempty(times)
    error('model and time should be specifiy')
end

h = waitbar(0,'Processing...','Name','Saving');    
for k = 1:length(times)
    waitbar(k/length(times),h)
    time = times{k};
    for m = 1:length(models)
        model = models{m};   
        LoadData

        %combined observation/regression matrices
        B_Lpj = blkdiag(B_LPJ,B_LPJ);
        A_Lpj = blkdiag(A_LPJ,A_LPJ);
        A1 = [B_Lpj, A_Lpj];

        Full_name = ['MCMC_',model,'_',num2str(time),'.mat'];
        filename = fullfile(pathname,Full_name);
        Full1 = load(filename);
        nu = mean(A1*[Full1.MCMC.beta(:,burnin:skip:end);Full1.MCMC.x(:,burnin:skip:end)],2);
        Full_MCMC = invalr( reshape(nu, [length(LPJ),d]) );
        results.(model).(time) = Full_MCMC;
        
        clearvars -except results times models time LPJ pathname skip burnin flag_save_csv h
    end
    if flag_save_csv
        Names_model = {'Lon','Lat','C_intercept','B_intercept','U_intercept','C_elev','B_elev','U_elev',...
                       'C_L_K_ESM','B_L_K_ESM','U_L_K_ESM','C_L_K_RCA3','B_L_K_RCA3','U_L_K_RCA3',...
                       'C_L_H_ESM','B_L_H_ESM','U_L_H_ESM','C_L_H_RCA3','B_L_H_RCA3','U_L_H_RCA3'};
        res = [LPJ(:,1:2),results.intercept.(time),results.elevation.(time),results.LPJ-KK10-ESM.(time),...
                          results.LPJ-KK10-RCA3.(time),results.LPJ-HYDE-ESM.(time),results.LPJ-HYDE-RCA3.(time)];

        name = ['results_',num2str(time),'.csv'];
        fullname = fullfile(pathname,name);
        fid = fopen(fullname,'w');
        fprintf(fid, [Names_model{1} sprintf(',%s',Names_model{2:end}) '\n']);
        fclose(fid);
        dlmwrite(fullname, res, '-append', ...
                 'delimiter', ',','precision','%.10f'); 
    end
    save(fullfile(pathname,'results_all.mat'),'results')
end
close(h)
end
