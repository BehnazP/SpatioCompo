function [ax,h]=suptitle(text)
% SUPTITLE centers a title over a group of subplots.
% Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.
%
% This code downloaded as an extension from MATLAB extenson community but
% unfortunately does not include writer information.

ax=axes('Units','Normal','Position',[.075 .1 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on');
title(text);
if (nargout < 2)
    return
end
h=get(ax,'Title');
end
