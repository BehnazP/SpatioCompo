if isunix
    addpath('./composition')
    addpath('./derivatives')
    addpath('./MCMC')
    addpath('./General')
    addpath('./Q')
    addpath('./Data')
elseif ispc
    addpath('.\composition')
    addpath('.\derivatives')
    addpath('.\MCMC')
    addpath('.\General')
    addpath('.\Q')
    addpath('.\Data')
end

%% run the mcmc
iter = 1e5;
delta = 1e-5;


times = {'AD1900','AD1725','BC4000'};
models = {'intercept','elevation','LPJ_KK10_ESM','LPJ_KK10_RCA3','LPJ_HYDE_ESM','LPJ_HYDE_RCA3'};

flag_save = 1;% save the results
mkdir Results%creat a folder to save the reuslts in

if ~exist('pathname','var')
    if isunix
        pathname = './Results';% the path for the results
        addpath(pathname)
    elseif ispc
        pathname = '.\Results';% the path for the results
        addpath(pathname)
    end
end


for m = 1:length(models)
    model = models{m};
    for j = 1:length(times)
        time = times{j};
        LoadData
        Main
    end
end

%% all the reconstructed

times = {'AD1900','AD1725','BC4000'};
models = {'intercept','elevation','LPJ_KK10_ESM','LPJ_KK10_RCA3','LPJ_HYDE_ESM','LPJ_HYDE_RCA3'};
burnin = 1e4;
skip = 10;

flag_save_csv = 0;% no results will be saved in csv format
SaveResults(times,models,pathname,burnin,skip,flag_save_csv)

%% plot results & EFI-FM

times = {'AD1900','AD1725','BC4000'};
LCC = {'Coniferous','Broadleaved','Unforested'};
titles = {'CE 1900', 'CE 1725', 'BCE 4000'} ;
models = {'intercept','elevation','LPJ_KK10_ESM','LPJ_KK10_RCA3','LPJ_HYDE_ESM','LPJ_HYDE_RCA3'};
load('results_all.mat')

flag_save_fig = 1; % 1 if you want to save the figure in pdf format
plotResults(results,times,titles,LCC,pathname,flag_save_fig);
