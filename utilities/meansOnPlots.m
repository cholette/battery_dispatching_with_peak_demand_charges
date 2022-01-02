function meansOnPlots(hf,means,stds)

N = length(hf.Children);
for ii = 1:N
    ax = hf.Children(ii);
    YY = ax.YLim;
    XX = ax.XLim;
    hold(ax,'on')
    M = means(ii);
    S = stds(ii);
    plot(ax,M*[1 1],YY,'r--','linewidth',2)
    str = {"Mean = "+num2str(M,'%.1f'),"STD = "+num2str(S,'%.1f')};
%     annotation('textbox',[0.8 0.8 0.1 0.1],'String',str,'FitBoxToText','on')
    text(ax,M+0.1*diff(XX),mean(YY),str,'HorizontalAlignment','left')
    hold(ax,'off')
end