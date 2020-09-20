% A script to generate the movie

close all
frame = 1;
if sim.displayCode==0,
    A_max = ceil(com.nSpecies);
else
    A_max = sim.displayAmax;
    Amax = 2;
end

for t=1:sim.displayInterval:sim.iter,

    A = squeeze(movieObj(frame,:,:));

    if t==1,
        % Create figure window

        hFig = figure(sim.figureNum);
        set(hFig,'Name',['SimCom ' date()],'Position',[10   383   514   559]);

        % Plot species richness
        if sim.displayCode==0 || sim.displayCode==-2,
            colormap(jet(ceil(A_max)+1))
            xTickLabel = round(linspace(0,A_max,6));
            xTick = 1 + xTickLabel;
        elseif sim.displayCode==-3,
            map = colormap(jet(com.nSpecies));
            hLine = plot([1:com.nSites(1)],full(A'),'LineWidth',2);

            for s=1:com.nSpecies,
                set(hLine(s),'Color',map(s,:))
            end
            hold on
            plot([1, com.nSites(1)],[1 1]*com.Amin,'k--','LineWidth',2)
            hold off

            axis([0 com.nSites(1)+1 0.9*com.Amin A_max]);

            set(gca,'FontSize',20)
            xlabel('Spatial Position')
            ylabel('Species abundance')
        else
            cmap = colormap(jet(128));
            cmap(1,:) = [1 1 1];
            colormap(cmap)
            xTickLabel = linspace(0,A_max,6);
            xTick = linspace(1,128,6);
        end

        if sim.displayCode~=-3,
            hImage = image( 1 + (length(get(gcf,'Colormap'))-1) * (A / A_max));
            set(gca,'XTick',[],'YTick',[]);
            % Add in a colour bar
            hCBar = colorbar('location','SouthOutside','XTick',xTick,'XTickLabel',xTickLabel);
            if sim.displayCode==0,
                set(get(hCBar,'XLabel'),'String','Species richness','FontSize',18)
            elseif sim.displayCode==-1
                set(get(hCBar,'XLabel'),'String','Mean abundance')
            elseif sim.displayCode==-2
                set(get(hCBar,'XLabel'),'String','Most abundant species id')
            elseif  sim.displayCode>0 & sim.displayCode<=com.nSpecies,
                set(get(hCBar,'XLabel'), ...
                    'String',['Abundance of species ' num2str(sim.displayCode)],...
                    'FontSize',18)
            end
        end
        title(['t = ' num2str(t)],'FontSize',18);
        
    else
        % Update figure
        if sim.displayCode~=-3,
            set(hImage,'cData',1 + (length(get(gcf,'Colormap'))-1) * (A/A_max))
        else
            for s=1:com.nSpecies,
                % Updata the data for each of the lines
                set(hLine(s),'YData',A(s,:))
            end
        end
        
        set(get(gca,'Title'),'String',['t = ' num2str(t)])
    end

    drawnow
    
    % Create the movie
    movieFinal(frame) = getframe(gcf);
    
    
    frame = frame + 1;
end
