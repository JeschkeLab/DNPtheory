function  fig=plot_IF_trajectory_3D( vec,rot_ax,Rtot,repeats,fig_num,mode,view_angle,frame_ind)


npoints = size(vec,2);
for rr = 1:repeats-1
    for ii = 1:npoints
        vec(:,ii+rr*npoints) = (Rtot')^rr*vec(:,ii);
    end
end


x=vec(1,:);
y=vec(2,:);
z=vec(3,:);


%calculate effective rot traj
n_eff_traj=100;
rotm_eff = rotaxi2mat(rot_ax,2*pi*320/360/n_eff_traj);
eff_traj_vec = [0; 0; 1];
for ii=2:n_eff_traj
   eff_traj_vec(:,ii)=rotm_eff*eff_traj_vec(:,ii-1); 
end



fig=figure(fig_num);
clf
hold on
% axis([-1 1 -1 1 -1 1])
coordsys_x=1.5*[0 1 -1 0 0 0 0 0 0 0];
coordsys_y=1.5*[0 0  0 0 1 -1 0 0 0 0 ];
coordsys_z=1.5*[0 0 0 0 0 0 0 1 -1 0];
plot3(coordsys_x,coordsys_y,coordsys_z,'k')
text(1.8,0,0,strcat('$S_x',frame_ind,'$'),'fontsize',15,'interpreter','latex')
text(0,1.8,0,strcat('$S_y',frame_ind,'$'),'fontsize',15,'interpreter','latex')
text(0,0.2,1.5,strcat('$S_z',frame_ind,'$'),'fontsize',15,'interpreter','latex')

% [x_sph y_sph z_sph] = sphere(50);
% r_sph = 0.99;
% h_sph = surfl(r_sph*x_sph, r_sph*y_sph, r_sph*z_sph);
% colormap gray
% set(h_sph, 'FaceAlpha', 0.2)
% shading interp

view([45 45])
set(gca,'Ztick',[],'Ytick',[],'Xtick',[],'Xcolor','w','Ycolor','w','Zcolor','w')
set(gcf, 'Color', 'w');
axis vis3d


%% Plot effective field vector if R is provided
if nargin >= 3
    if  rot_ax'*rot_ax>0.5
        stemWidth = 0.025;
        if rot_ax(3)>0
            mArrow3([0 0 0],2.5*rot_ax,'color',[0.8 0.2 0.2],'stemWidth',stemWidth);
        else
            mArrow3([0 0 0],-2.5*rot_ax,'color',[0.8 0.2 0.2],'stemWidth',stemWidth);
        end
    end
end


%% static plot
if(strcmp('static',mode))
    
    plot3(eff_traj_vec(1,:),eff_traj_vec(2,:),eff_traj_vec(3,:),'color',0.5*[1 1 1])
    mArrow3(eff_traj_vec(:,end-3),eff_traj_vec(:,end),'color',0.5*[1 1 1]);
    
    plot3(x(npoints:end),y(npoints:end),z(npoints:end),'b','linewidth',2)
    plot3(x(1:npoints),y(1:npoints),z(1:npoints),'r','linewidth',3)
    for ii=1:repeats
        plot3(x((ii)*npoints),y((ii)*npoints),z((ii)*npoints),'ok','markerfacecolor','k')
    end
end



%% dynamic plot
figure(fig);
drawnow
if(strcmp('film',mode))
    h=plot3(x(1:1),y(1:1),z(1:1));
    
    for ii=2:length(x)
        
        %only plot if something actually changes, otherwise one mainly
        %waits
        if(x(ii)==x(ii-1) && y(ii)==y(ii-1) && z(ii)==z(ii-1))
        else
            %         delete(h)
            h=plot3(x(1:ii),y(1:ii),z(1:ii),'b');
            pause(0.001)
        end
        
    end
    
end


view(view_angle)
fig.Renderer = 'painters';

axis tight

end

