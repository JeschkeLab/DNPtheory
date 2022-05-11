function h=IF_plot(time,R,fig_num)


%% Plot IF trajectories of 
% h=figure('Units','pixels','Position',[114.6000 153.8000 1200 400.2000]);
h=figure(fig_num);
set(h,'units','points')
set(h,'Position',[250 150 420 360])

clf
subplot(3,3,1)
IF_plot_ij(time,squeeze(R(1,1,:)),'r');
text(time(end)/2,1.4,'$x$','horizontalalignment','center','interpreter','latex')
ylab=ylabel('$\widetilde{x}$','rotation',0,'horizontalalignment','right','verticalalignment','middle','fontsize',12,'interpreter','latex');
ylab.Position=ylab.Position;
subplot(3,3,4)
IF_plot_ij(time,squeeze(R(2,1,:)),'r');
ylab=ylabel('$\widetilde{y}$','rotation',0,'horizontalalignment','right','verticalalignment','middle','fontsize',12,'interpreter','latex');
ylab.Position=ylab.Position;subplot(3,3,7)
IF_plot_ij(time,squeeze(R(3,1,:)),'r');
ylab=ylabel('$\widetilde{z}$','rotation',0,'horizontalalignment','right','verticalalignment','middle','fontsize',12,'interpreter','latex');
ylab.Position=ylab.Position;

subplot(3,3,2)
IF_plot_ij(time,squeeze(R(1,2,:)),[0 0.7 0]);
text(time(end)/2,1.4,'$y$','horizontalalignment','center','interpreter','latex')
subplot(3,3,5)
IF_plot_ij(time,squeeze(R(2,2,:)),[0 0.7 0]);
subplot(3,3,8)
IF_plot_ij(time,squeeze(R(3,2,:)),[0 0.7 0]);


subplot(3,3,3)
IF_plot_ij(time,squeeze(R(1,3,:)),'b');
text(time(end)/2,1.4,'$z$','horizontalalignment','center','interpreter','latex')
subplot(3,3,6)
IF_plot_ij(time,squeeze(R(2,3,:)),'b');
subplot(3,3,9)
IF_plot_ij(time,squeeze(R(3,3,:)),'b');

for ii=1:numel(h.Children)
    h.Children(ii).XTick=[];
    h.Children(ii).YTick=[];
end


set(findall(h,'-property','LineWidth'),'LineWidth',1.5); % change linewidth


end

function h=IF_plot_ij(x,y,color)
mygray = 0.7*[1 1 1];

hold on
axis([0 x(end) -1 1])
plot(xlim,[0 0],'color',mygray)
plot(x(end)/4*[1 1],ylim,'color',mygray)
plot(x(end)/2*[1 1],ylim,'color',mygray)
plot(x(end)*3/4*[1 1],ylim,'color',mygray)
plot(x,y,'color',color)
plot(x([1 end]),y([1 end]),'o','color',color,'markerfacecolor',color)
box on

end