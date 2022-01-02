function visualizeTraining(k,tk,Xsk,phis,Vk,rbfVFA)

fh = findobj( 'Type', 'Figure', 'Name', "Estimation" );
if isempty(fh)
    figure("Name","Estimation")
else
    figure(fh)
end

if iscell(rbfVFA.xi)
    xi = rbfVFA.xi{k};
    bar(phis*rbfVFA.theta{k})
else
    xi = rbfVFA.xi;
    bar(phis*rbfVFA.theta(:,k))
end

hold on
plot(Vk,'*')
drawnow
hold off

fh = findobj( 'Type', 'Figure', 'Name', "Value Function" );
if isempty(fh)
    figure("Name","Value Function")
else
    figure(fh)
end
scatter3(Xsk(1,:),Xsk(2,:),Xsk(3,:),[],Vk)

xlim([min(xi(1,:)),max(xi(1,:))])
ylim([min(xi(2,:)),max(xi(2,:))])
zlim([min(xi(3,:)),max(xi(3,:))])
xlabel('demand')
ylabel('soc')
zlabel('recorded peak')
title("k="+num2str(k)+", tod="+num2str(hour(tk),'%02d')+":"+num2str(minute(tk),'%02d'));
colorbar
drawnow
hold off