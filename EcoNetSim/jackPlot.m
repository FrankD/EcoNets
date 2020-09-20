% Plot food web topology
close all


figure(1)
set(gcf,'Name','Food web topology')

levels = unique(model.trophicLevel);

colormap(pink(1+numel(levels)));
web = repmat(full(model.trophicLevel),1,com.nSpecies)';
web(find(model.webTopology)) = -1;
imagesc(web)
axis square
xlabel('Species ID (consummers)')
ylabel('Species ID (prey)')
set(gca,'XAxisLocation','top')
hc1 = colorbar('Location','SouthOutside');

stepSize = (max(levels)+1)/(numel(levels)+1);
set(hc1,'XLim',[-1 + stepSize,max(levels)])
ticks = -1+stepSize + stepSize*[0.5:1:numel(levels)];
set(hc1,'XTick',-0.25 + 0.75*[0.5,1.5,2.5])
labels = num2cell(levels);
labels{1} = 'Basal';
set(hc1,'XTickLabel',labels)
set(get(hc1,'XLabel'),'String','Minimum Trophic Level')

figure(2)
set(gcf,'Name','Species richness')
colormap(hot(range(richness(1:end))+1))

imagesc(richness)
axis square
title('Species richness')
hc2 = colorbar('Location','EastOutside');
%set(get(hc2,'YLabel'),'String','Species richness')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
set(gcf,'Name','Most abundant species')

colormap(jet(com.nSpecies))
imagesc(domID)
axis square
set(gca,'XAxisLocation','bottom')
xlabel('X coordinate')
xlabel('Y coordinate')
title('Most abundant species')

%hc3 = colorbar('Location','SouthOutside');
%set(get(hc3,'XLabel'),'String','Most abundant species ID')
hold on
plot(round(com.nSites(1)/2)*[1 1],[2 com.nSites(2)-1],'k--','LineWidth',2)
text(round(com.nSites(1)/2)*[1 1],[1 com.nSites(2)],{'A','A'}, ...
    'HorizontalAlignment','Center','VerticalAlignment','Middle')
plot([2 com.nSites(1)-1],round(com.nSites(2)/2)*[1 1],'k--','LineWidth',2)
text([1 com.nSites(1)],round(com.nSites(2)/2)*[1 1],{'B','B'}, ...
    'HorizontalAlignment','Center','VerticalAlignment','Middle')
hold off
map = colormap();


figure(4)
set(gcf,'Name','Abundance transects')

ind = false(com.nSites);
ind(:,round(com.nSites(2)/2)) = true;
transect1 = A(:,reshape(ind,1,prod(com.nSites)));
ind = false(com.nSites);
ind(round(com.nSites(1)/2),:) = true;
transect2 = A(:,reshape(ind,1,prod(com.nSites)));

subplot(2,1,1)
for i=1:com.nSpecies,
    h1(i) = plot([1:com.nSites(2)],transect1(i,:));
    hold on
    set(h1(i),'Color',map(i,:))
end
axis([1 com.nSites(2) 0 max(transect1(1:end))+0.5])
ticks = linspace(1,com.nSites(2),5);
ticks(2:end-1) = 5*round(ticks(2:end-1)/5);
label = num2cell(ticks);
label{1} = 'A';
label{end} = 'A';
set(gca,'XTick',ticks,'XTickLabel',label);


subplot(2,1,2)
for i=1:com.nSpecies,
    h2(i) = plot([1:com.nSites(1)],transect2(i,:));
    hold on
    set(h2(i),'Color',map(i,:))
end
axis([1 com.nSites(1) 0 max(transect2(1:end))+0.5])
ticks = linspace(1,com.nSites(1),5);
ticks(2:end-1) = 5*round(ticks(2:end-1)/5);
label = num2cell(ticks);
label{1} = 'B';
label{end} = 'B';
set(gca,'XTick',ticks,'XTickLabel',label);

xlabel('Spatial Location')
ylabel('Species abundance')
