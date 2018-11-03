function save_fig(SaveName,orientation)
% SAVE_FIG take the current figure handle and save the figure in .pdf format
% 
% save_fig(SaveName,orientation)
% input:
% SaveName          name of the file
% orientation       "portrait" or "landscape"
%
% SAVE_FIG.m 2018-07-15 Behnaz@pirzamanbin.name$

h=gcf;
set(h,'PaperOrientation',orientation);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', SaveName);  

end