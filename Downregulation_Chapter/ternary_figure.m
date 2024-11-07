
% advanced_example.m This script shows a more involved case of
% ternary_plots, where the three axes are not space 0->1 & do not sum to 1.
% Additional options are added to show how various elements can be
% customized. 
% clear all; close all; clc

function ternary_figure(A,B,Z,wlimits,custom_labels,climit)
arguments
    A 
    B 
    Z 
    wlimits 
    custom_labels
    climit = [0 1]
end
showMaxMin=false;
changeAxesColor=true;

ccmap=cmocean('deep',80);
cmap=ccmap(1:end-10,:);
%% (1) Add paths 
    % This function can be run automatically if copied into "userpath/sartup.m
    add_ternary_paths
    
%% (2) Create Figure with two subplots
    % ffig = figure('Name','Advanced Example','Position',[100 100 1000 400]);
   
     % wlimits = ternary_axes_limits( 100,'l',20,'low',...
     %                                   'l',80,'high',...
     %                                   'r',10,'low', false ); % turn on example plot

     % wlimits = ternary_axes_limits(1,'l',0,'low','l',1,'high','r',0,'low', false);

%% (5) Create the Axes using customized settings; 
    
    % "vgen" is a cell array of custom settings specific to Ternary_Plots.
    %    All options are listed in ternary_axis.m ->  initialize_ternary_handle(). 
    %    This is an example of extreme customization.
    vgen  = { 'wlimits',       wlimits ,... % Axes will match wlimits ranges
              'gridspaceunit', 5,      ... % Number of grid lines
              'ticklinelength', [.01,.01,.01],    ... % [3.25,3.25,4]length of tick-lines in ABC coord.
              'tick_fmt',      '%.1f', ... % tick label formatting
              'titlelabels', custom_labels, ... % custom labels custom_labels={'j_N','j_L','j_{Si}'}
              'titlerotation', [0,0,0], ... % Set all titles to horizontal
              'link_color', {'tick','title'},... % Link all axes colors
              'titleshift',[ -0.25, 0, 0.25; 0.09, -0.18, 0.09 ]... % shift titles [ -0.18, 0, 0.18; 0.085, -0.11, 0.085 ]
              'tickshift', [-0.02, -0.02, 0.0; -0.01,-0.00,0]
            };
    
	% Ternary Axes Outline   - Passed directly to plot3() 
    vout  = { 'LineWidth', 2, 'LineStyle', '-','Color','k'};
        
    % Ternary gridlines  - Passed directly to plot3()
    vgrid = { 'LineStyle','-','LineWidth',1, 'Color',[0 0 0 0.2] };
        
    % Ternary Tick Line - Passed directly to plot3()
    vtick_line = { 'LineStyle', '-', 'LineWidth', .5, 'Color',[0 0 0 1.0] };
    
    % Ternary Tick Labels - Passed directly to text()
    vtick_label = { 'FontWeight','normal', 'FontSize', 8 };
       
    % Ternary Axes Label  - Passed directly to text()
    vlab  = { 'FontWeight','normal', 'FontSize', 10 };
    
    % First time a plot appears
    % Create Ternary Axes & return Handle
    handle = ternary_axes( vgen, vout, vgrid,vtick_line, vtick_label, vlab );   

%% (8) Plot the surface Z defiend on rows of ABC points in xmat
    
    % Define a customized color bar
    Cbar = {'FontWeight','normal'};%,'Position',[0.90 0.17 0.02 0.68] };
    
    % Create Surface Plot + colorbar
    dataplots(1).obj = ternary_surf( wlimits, 'l', A, 'b', B, Z , Cbar );
    
    % Set shading (e.g. flat or interp)
    shading('interp')
    
    % Set Colormap to B/W
    colormap(cmap)
    
    % Reset Range of Surface Colors
    clim(climit)
    
    % Add custom Data tip
    set(datacursormode(gcf),'UpdateFcn',{@ternary_datatip,dataplots(1).obj.ZData,wlimits})
    
    % Add Title
    % title('A Fruity Example','FontSize',18)
    
    if showMaxMin==true
        % (9) Add labels at Max
        [val,idx] = max(Z(:));
        str = ['Max=',num2str(val,'%4.1f')];
        var = {'MarkerFaceColor','b','MarkerEdgeColor','b'};
        dataplots(2).obj = ternary_scatter3(wlimits, 1, A(idx), 2, B(idx), [], 'none', var{:});
        var = {'Color','b'};
        dataplots(3).obj = ternary_text(wlimits, 1, A(idx), 2, B(idx),str, [] , var{:} );
        dataplots(3).obj.HorizontalAlignment ='center';
        dataplots(3).obj.VerticalAlignment   ='top';

        %(10) Add labels at Min
        [val,idx] = min(Z(:));
        str = ['Min=',num2str(val,'%4.1f')];
        var = {'MarkerFaceColor','c','MarkerEdgeColor','c'};
        dataplots(4).obj = ternary_scatter3(wlimits, 1, A(idx), 2, B(idx), [], 'none', var{:});
        var = {'Color','c'};
        dataplots(5).obj = ternary_text(wlimits, 1, A(idx), 2, B(idx),str, [] , var{:} );
        dataplots(5).obj.HorizontalAlignment ='center';
        dataplots(5).obj.VerticalAlignment   ='bottom';
    end
%% (9) Adjust Colors of axes
if changeAxesColor==true
    % Change the Apples title, tick label, tick lines color
    % handle.link_color = {'tick','title'};
    % handle = adjust_axis_color(handle,'left',[0.6353, 0.0784, 0.1843, 0.6]);

    % Change LineStyle of grid lines
    handle.grid.lines(1,1).LineWidth = 1.5;
    handle.grid.lines(1,1).LineStyle = ':';
    handle.grid.lines(1,2).LineWidth = 1.0;
    handle.grid.lines(1,3).LineStyle = '--';

    % Change color of Oranges grid lines and tick labels
    handle.link_color = {'grid','tick_label'};
    handle = adjust_axis_color(handle,'bot',[0.85, 0.33, 0.098, 0.7 ] );

    % Change color of Banana title
    % handle.link_color = {'title'};
    % handle = adjust_axis_color(handle,'r',[0.58, 0.54, 0.0] );

    % Restack Final Plots
    handle = restack_dataplots(handle, dataplots);
end
    
    