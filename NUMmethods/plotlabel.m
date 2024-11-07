%
% Add a label to a plot in the top left corner
% if bBackground = true then the background is blotted out.
%
function h = plotlabel(sLabel,bBackground)
if (nargin==1)
    bBackground = true;
end

heightlabel = 0.18; % Height of label in cm

% Get the size of the paper:
set(gcf,'paperunits','centimeters')
set(gcf,'units','centimeters')
tmp = get(gcf,'paperposition');
widthpaper = tmp(3);
heightpaper = tmp(4);

% Get relative size of the plot:
pos = get(gca,'position');
height = pos(4);
width = pos(3);

x = heightlabel/(widthpaper*width);
y = (heightpaper*height-0.2*heightlabel)/(heightpaper*height);

%text(x,y,sLabel,'fontsize',10,'fontname','arial','fontweight','bold','units','normalized','color',[1 1 1])
h=text(x,y,sLabel,'fontsize',10,'fontname','arial',...
  'units','normalized',...
  'verticalalignment','top');
if bBackground
    set(h,'backgroundcolor',[1 1 1])
end
