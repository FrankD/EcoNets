function [A, model, movieObj] = communitySimulation_2d(sim,com),
% This function runs the main community simulation program

% This version uses a 2d spatial arena
% The spatial data is rearranged so that it it one long vector

%%% Generate the basic model variables
model = generateModel_2d(sim,com);
A = full(model.initialA);

%%% The main loop
if sim.movie,
    frame = 1;
    nFrames = 1 + (sim.iter-1-mod(sim.iter-1,sim.displayInterval)) / ...
        sim.displayInterval;
    movieObj = zeros(nFrames,com.nSites(1),com.nSites(2));
else
    movieObj = 0;
end

locations = [];


for i=1:sim.iter,  % Loop through all the iterations
    if mod(i-1,sim.displayInterval)==0,
        if sim.debug
          disp(['Iteration number = ' num2str(i)])
        end
        if com.nSites>1 & sim.display
            imageData = displayData(i,A,com,sim);

            if sim.movie,
                movieObj(frame,:,:) = imageData;
                frame = frame+1;
            end

        end
    end

    if sim.changeGrowthRate & i==5,
        disp('Changing growth rates')
        Y = reshape(repmat([1:com.nSites(1)]'-(com.nSites(1)+1)/2,com.nSites(1),1),1,prod(com.nSites));
        for s=1:com.nSpecies,
            model.speciesGrowthRates(s,:) = model.speciesGrowthRates(s,:) + Y;
        end
    end

    if sim.changeInteractions,
        if i<10,
            weight = 0;
        else
            weight = min([(i-10)/40,1]);
        end
    else
        weight = 1;
    end

    r = full(model.speciesGrowthRates);
    gam = full(model.speciesIntraDD);
    
    A_disp = A * model.dispersalMatrix;
    
    demoStoch = sqrt(com.e_sigma)*randn(size(A_disp));    
    
    intraDD = intraDensityDependence(A_disp, gam);
    
    interactions = weight * model.speciesInteractionStrengths;
    
    interDD = interDensityDependence(A_disp, interactions);
    
    A =  A_disp .* (full(spones(r)) + ...
            (r - interDD - intraDD + demoStoch)*sim.dt);

    extinctThreshold = -com.Amin * log(rand(size(A)));
    
    A(A - extinctThreshold <= 0) = 0;
    
    if(i == 100)
      [same, total, locations] = count_outliers(A, locations);
    end
    
    if(mod(i, 100) == 0)
      %real_interactions = full(spones(model.speciesInteractionStrengths));
      %display(score_of_scores(standardizeCols(A')', real_interactions));
      %[same, total, new_locations] = count_outliers(A, locations);
      %display(same);
      %display(total);
%      imagesc(locations);
 %     pause(5);
      
    end
    
%     for n=1:prod(com.nSites),  % Loop through each resource patch
%         r_Site = full(model.speciesGrowthRates(:,n));
%         gam_Site = full(model.speciesIntraDD(:,n));
%         
%         if mod(n, round(prod(com.nSites)/20)) == 0 && sim.debug
%           fprintf('Update done %g\n', n/prod(com.nSites))
%         end
%         
%         % Calculate species abundance after dispersal for the nth site
%         % model.dispersalMatrix(i,n) is the prob that an individual from
%         % site i ended up at site n
%         %        tmp = real(A) * spdiags(model.dispersalMatrix(:,n),0,prod(com.nSites),prod(com.nSites));
%         An_Site = A * model.dispersalMatrix(:,n);
%         %        tmp(tmp<1e-10) = 0;
% 
%         % Calculate density of dispersing species
%         %        An_disperse = sparse(tmp);
% 
%         % Calculate new abundances using Engen and Lande's model
%         %        An_Site = sum(An_disperse,2);
%         demoStoch = sqrt(com.e_sigma)*randn(size(An_Site));
% 
%         % Calculate inter-species density dependence
%         intraDD = intraDensityDependence(An_Site,gam_Site);
% 
% 
%         interactions = weight * model.speciesInteractionStrengths;
%         interDD = interDensityDependence(An_Site,interactions);
% 
%         An_Site =  An_Site .* (full(spones(r_Site)) + ...
%             (r_Site - interDD - intraDD + demoStoch)*sim.dt);
% 
%         % Keep those species who haven't gone extinct
%         extinctThreshold = -com.Amin * log(rand(size(An_Site)));
%         %An_Site(An_Site>eps & An_Site-extinctThreshold<=0) = 0;
%         An_Site(An_Site-extinctThreshold<=0) = 0;
%         % Calculate new abundance
%         A_new(:,n) = sum(An_Site,2);
%     end
% 
%     % Replace old abundance with a new one
%     A = A_new;

    % Calculate number of species
    num_species(i,:) = sum(A>10^5*eps);
    tot_num_species(i) = sum(sum(A,2)>10^5*eps)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [same, number, locations] = count_outliers(A, previous) 
  mu = mean(A, 2);

  A_temp = A - repmat(mu, 1, size(A, 2));

  stand = std(A, 0, 2);
  counter = 0;
  
  locations = zeros(size(A));
  
  for i = 1:size(A, 1)
    A_row = A_temp(i, :);
  
    outliers = find(abs(A_row) > 2.75*stand(i));  
    
    locations(i, outliers) = 1;

    counter = counter + size(outliers, 2);
  end

  
if(size(previous) == [0, 0])
  same = counter;
else
  previous(previous == 0) = -1;
  same = sum(sum(previous == locations));
end
  
number = counter;



%% Interspecific Density dependence
function interDD = interDensityDependence(Aj_Site,I),
% Calculate inter-species density dependence
interDD = I * Aj_Site;
%interDD(interDD<1e-10) = 0;
return



%% Intraspecific Density dependence
function intraDD = intraDensityDependence(Aj_Site,gam),
% Calculate intra-species density dependence

% Ths is a gompertz form
intraDD = gam .* spfun(@log, Aj_Site);
return

%% Display data
function A = displayData(t,X,com,sim),
% Function that displays the species abundance data


if sim.displayCode ==0,
    A = full(reshape(sum(X>com.Amin,1),com.nSites));
    A_max = com.nSpecies;
elseif sim.displayCode==-1
    A = full(reshape(mean(X,1),com.nSites));
    A_max = sim.displayAmax;
elseif sim.displayCode==-2
    [dum,A] = max(X,[],1);
    A = reshape(A,com.nSites);
    A_max = com.nSpecies;
elseif sim.displayCode==-3
    ind = false(com.nSites);
    ind(:,round(com.nSites(2)/2)) = true;
    A = X(:,reshape(ind,1,prod(com.nSites)));
    A_max = sim.displayAmax;
elseif sim.displayCode>0 & sim.displayCode<=com.nSpecies,
    A = full(reshape(X(sim.displayCode,:),com.nSites));
    A_max = sim.displayAmax;
end

if ~exist('hFig'),
    % Create figure window
    persistent hFig hLine
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
        %        hLine = get(hf,'Children');
        %        hLine = findobj('Type','Line');
        for s=1:com.nSpecies,
            set(hLine(s),'Color',map(s,:))
        end
        hold on
        plot([1, com.nSites(1)],[1 1]*com.Amin,'k--','LineWidth',2)
        hold off
        %        y_max = 10^full(log10(ceil(max(nonzeros(A)))));
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
        %        set(gca,'XAxisLocation','top')
        set(gca,'XTick',[],'YTick',[]);
        % Add in a colour bar
        hCBar = colorbar('location','SouthOutside','XTick',xTick,'XTickLabel',xTickLabel);
        if sim.displayCode==0,
            set(get(hCBar,'XLabel'),'String','Species richness')
        elseif sim.displayCode==-1
            set(get(hCBar,'XLabel'),'String','Mean abundance')
        elseif sim.displayCode==-2
            set(get(hCBar,'XLabel'),'String','Most abundant species id')
        elseif  sim.displayCode>0 & sim.displayCode<=com.nSpecies,
            set(get(hCBar,'XLabel'),'String',['Abundance of species ' num2str(sim.displayCode)])
        end
    end
    title(['t = ' num2str(t)]);

else
    % Update figure
    if sim.displayCode~=-3,
        set(hImage,'cData',1 + (length(get(gcf,'Colormap'))-1) * (A/A_max))
    else
        for s=1:com.nSpecies,
            % Updata the data for each of the lines
            set(hLine(s),'YData',full(X(s,:)))
        end
    end

    set(get(gca,'Title'),'String',['t = ' num2str(t)])
end

drawnow
return
