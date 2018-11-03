function plotResults(results,times,titles,LCC,pathname,flag_save)
% PLOTRESUTS plot the reconstructed map from MCMC samples and have an
%            option to save them in .pdf
%
% plotResults(results,times,titles,LCC,pathname,flag_save)
% input:
% results       a structure containing reconsturcted maps from MCMC sample
% times         different time periods
% titles        different titles for the suptitle
% LCC           land cover compositions names
% pathname      path for saving plots
% flag_save     indicator for saving the plots, if 1 the plot will be saved
%
% PLOTRESULTS.m 2018-07-15 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

model = 'elevation';
h = waitbar(0,'Processing...','Name','Saving');

for t = 1:length(times)
    waitbar(t/length(times),h)
    time = times{t};
    ti = titles{t};
    LoadData

    clearvars -except LPJ results Data times t time A_LPJ A_REV sz REV Rev LCC titles ti I_LPJ I_Rev model pathname flag_save h

    REV = A_REV'*Rev(:,3:5);

    I_LPJ = reshape(A_LPJ'*I_LPJ',[sz,1]);
    I_Rev = reshape(A_REV'*I_Rev',[sz,1]);

    [EFMadj, ~]=loc(LPJ(:,1:2),Data.EFMadj1);
    EFI = reshape(A_LPJ'*EFMadj(:,3:5),[sz 3]);

    res_intercept = reshape(A_LPJ'*results.intercept.(time),[sz,3]);
    res_elevation = reshape(A_LPJ'*results.elevation.(time),[sz,3]);
    res_LPJ_KK10_ESM = reshape(A_LPJ'*results.LPJ_KK10_ESM.(time),[sz,3]);
    res_LPJ_KK10_RCA3 = reshape(A_LPJ'*results.LPJ_KK10_RCA3.(time),[sz,3]);
    res_LPJ_HYDE_ESM = reshape(A_LPJ'*results.LPJ_HYDE_ESM.(time),[sz,3]);
    res_LPJ_HYDE_RCA3 = reshape(A_LPJ'*results.LPJ_HYDE_RCA3.(time),[sz,3]);

    Rev = reshape(REV,[sz,3]);

    figure
    for k = 1:3
        if t == 1
            subplot(8,3,k)
            imagesc(Rev(:,:,k),'AlphaData',I_Rev)
            title(LCC(k))
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'PbLCC','rotation',90)
            end
            subplot(8,3,k+3)
            imagesc(EFI(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'EFI-FM','rotation',90)
            end
            subplot(8,3,k+6)
            imagesc(res_intercept(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'Intercept','rotation',90)
            end
            subplot(8,3,k+9)
            imagesc(res_elevation(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'Elevation','rotation',90)
            end
            subplot(8,3,k+12)
            imagesc(res_LPJ_KK10_ESM(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-KK10_{ESM}','rotation',90)
            end
            subplot(8,3,k+15)
            imagesc(res_LPJ_KK10_RCA3(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-KK10_{RCA3}','rotation',90)
            end
            subplot(8,3,k+18)
            imagesc(res_LPJ_HYDE_ESM(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-HYDE_{ESM}','rotation',90)
            end
            subplot(8,3,k+21)
            imagesc(res_LPJ_HYDE_RCA3(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-HYDE_{RCA3}','rotation',90)
            end
            h1 = get(subplot(8,3,24),'Position');
            colorbar('Position', [h1(1)+h1(3)+0.02  h1(2)  0.02  h1(2)+h1(3)*3.1])
        else

            subplot(7,3,k)
            imagesc(Rev(:,:,k),'AlphaData',I_Rev)
            caxis([0 1])
            title(LCC(k))
            axis xy
            axis off
            if k ==1
                text(-5,10,'PbLCC','rotation',90)
            end
            subplot(7,3,k+3)
            imagesc(res_intercept(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'Intercept','rotation',90)
            end
            subplot(7,3,k+6)
            imagesc(res_elevation(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'Elevation','rotation',90)
            end
            subplot(7,3,k+9)
            imagesc(res_LPJ_KK10_ESM(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-KK10_{ESM}','rotation',90)
            end
            subplot(7,3,k+12)
            imagesc(res_LPJ_KK10_RCA3(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-KK10_{RCA3}','rotation',90)
            end
            subplot(7,3,k+15)
            imagesc(res_LPJ_HYDE_ESM(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-HYDE_{ESM}','rotation',90)
            end
            subplot(7,3,k+18)
            imagesc(res_LPJ_HYDE_RCA3(:,:,k),'AlphaData',I_LPJ)
            caxis([0 1])
            axis xy
            axis off
            if k ==1
                text(-5,10,'LPJ-HYDE_{RCA3}','rotation',90)
            end
            h1 = get(subplot(7,3,21),'Position');
            colorbar('Position', [h1(1)+h1(3)+0.02  h1(2)  0.02  h1(2)+h1(3)*3.2])
        end
    end
    suptitle(ti);

    if flag_save
        SaveName = ['res_',num2str(time),'.pdf'];
        savename = fullfile(pathname,SaveName);
        save_fig(savename,'portrait')
    end
end
close(h)
end
