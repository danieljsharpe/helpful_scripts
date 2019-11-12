% generic script for a line plot from a .dat file in Matlab
clc, clf;
% absolute path to data files 1 and 2
abs_data_path = '/home/djs244/GROUP_DATABASES/LJ38/interbasin_cum';
fname1 = "rate_lim_cut1000.dat";
fname2 = "ngt_cum.dat";
ffmt1="%d %f %f %d";
ffmt2="%d %f %f %f %f";
col1_idx1 = 4; % column no. containing x data, file 1
col2_idx1 = 3; % column no. containing y data, file 1
col1_idx2 = 1; % column no. containing x data, file 2
col2_idx2 = 3; % column no. containing y data, file 2

% set dimensions of figure and of page
fig_height=400; fig_width=550; % in px
page_height=10; page_width=12.5; % in cm

fid=fopen(fullfile(abs_data_path,fname1));
celldata1=textscan(fid,ffmt1); % read data into cell array
fclose(fid);
xdatacell1=celldata1(:,col1_idx1);
ydatacell1=celldata1(:,col2_idx1);
xdata1=cell2mat(xdatacell1);
ydata1=cell2mat(ydatacell1);

fid=fopen(fullfile(abs_data_path,fname2));
celldata2=textscan(fid,ffmt2); % read data into cell array
fclose(fid);
xdatacell2=celldata2(:,col1_idx2);
ydatacell2=celldata2(:,col2_idx2);
xdata2=cell2mat(xdatacell2);
ydata2=cell2mat(ydatacell2);

yyaxis left % set active axis to left axis
plot(xdata1,ydata1,"r","Linewidth",3.);
hold on;
xlabel('Path number','Interpreter','latex','FontSize',14)
ylabel('Path cost','Interpreter','latex','FontSize',14)
%set the x limits and tick intervals here
xlimits=[0,1100];
ylimits1=[40.,80.];
ylimits2=[0.,2.5];
exponent=1.E27; % set to rescale data
ydata2=ydata2*exponent;
xtickintvl=200;
ytickintvl1=5;
ytickintvl2=0.5;
xlim([xlimits(1) xlimits(2)])
ylim([ylimits1(1) ylimits1(2)])
xticks(xlimits(1):xtickintvl:xlimits(2));
yticks(ylimits1(1):ytickintvl1:ylimits1(2));
yyaxis right
plot(xdata2,ydata2,"b","Linewidth",3.);
ylabel('Accumulated $k_{AB}^\mathrm{NGT}$ $\times~10^{27}$','Interpreter','latex','FontSize',14)
ylim([ylimits2(1) ylimits2(2)])
yticks(ylimits2(1):ytickintvl2:ylimits2(2));
lgd = legend('Path cost','Accumulated $k_{AB}^\mathrm{NGT}$','Interpreter','latex','FontSize',8,'location','southeast');
legend('boxoff');
ax = gca;
ax.YAxis(2).Exponent = 0;
set(ax,'TickDir','out');
set(ax,'TickLabelInterpreter','latex','FontSize',12)
set(ax,'YColor','k');
ytickformat('%.1f')
set(gcf,'position',[0,0,fig_width,fig_height]);
set(gcf,'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',[0,0,page_width,page_height]);
yyaxis left
ax = gca;
set(ax,'YColor','k');
set(ax,'TickLabelInterpreter','latex','FontSize',12)
% matlab saving syntax
saveas(gcf,fullfile(abs_data_path,'lj38_fi_pathcosts_ngtblock'),'epsc')
saveas(gcf,fullfile(abs_data_path,'lj38_fi_pathcosts_ngtblock'),'png')