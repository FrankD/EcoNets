function [X, deltaX, model, movieObj] = communitySimulation3(sim,com),
% This function runs the main community simulation program

%%% Generate the basic model variables
model = generateModel(sim,com);
X = model.initialX;

%%% The main loop
if sim.movie,
    frame = 1;
else
    movieObj = '';
end
    
for i=1:sim.iter,  % Loop through all the iterations
    if mod(i-1,sim.displayInterval)==0,
        disp(['Iteration number = ' num2str(i) ' Number of species ' num2str(nnz(any(X,2)))])
        if com.nSites>1 & sim.display
            displayData(i,X,com,sim);
            
                if sim.movie & i>1,
                    movieObj(frame) = getframe(gcf);
                    frame = frame+1;
                end
             
        end
    end

    for n=1:com.nSites,  % Loop through each resource patch
        r_Site = model.speciesGrowthRates(:,n);
        gam_Site = model.speciesIntraDD(:,n);
        % Calculate species abundance after dispersal for the nth site
        tmp = X * spdiags(model.dispersalMatrix(:,n),0,com.nSites,com.nSites);
        % Calculate density of dispersing species
        Xn_disperse = sparse(tmp);
        % Calculate new abundances using Engen and Lande's model
        Xn_Site = sum(Xn_disperse,2);
        demoStoch = sqrt(com.e_sigma)*sprandn(Xn_Site);

        % Calculate inter-species density dependence
        intraDD = intraDensityDependence(Xn_Site,gam_Site);
        interDD = interDensityDependence(Xn_Site/com.Xmin,model.speciesInteractionStrengths);

        Xn_Site =  Xn_Site .* (spones(r_Site) + ...
            (r_Site + interDD - intraDD + demoStoch)*sim.dt);

        % Keep those species who haven't gone extinct
        extinctThreshold = -com.Xmin * log(rand(size(Xn_Site)));
        Xn_Site(Xn_Site>eps & Xn_Site-extinctThreshold<=0) = 0;

        % Calculate new abundance and competitive abilities
        X_new(:,n) = sum(Xn_Site,2);
    end
    % Replace old abundance with a new one
    deltaX = sparse(X_new - X);
    % If species are independent then remove any species where groth rate
    % is negative across all space
 %   ind = find(max(deltaX,[],2)<-eps);
 %   X_new(ind,:) = 0;
    
    X = sparse(real(X_new));

    % Calculate number of species
    num_species(i,:) = sum(X>10^5*eps);
    tot_num_species(i) = sum(sum(X,2)>10^5*eps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Interspecific Density dependence
function interDD = interDensityDependence(Xj_Site,I),
% Calculate inter-species density dependence
interDD = I * Xj_Site;
%tmp = Xj_Site + spones(Xj_Site);

%tmp = Xj_Site;
%interDD = I * spfun(@log,tmp);
interDD(interDD<1e-10) = 0;
interDD = sparse(interDD);
return



%% Intraspecific Density dependence
function intraDD = intraDensityDependence(Xj_Site,gam),
% Calculate intra-species density dependence

% Ths is a gompertz form
%tmp = Xj_Site + spones(Xj_Site);
tmp = Xj_Site;
intraDD = gam .* spfun(@log,tmp);
return

%% Display data
function movieObj = displayData(t,X,com,sim,movieObj),
% Function that displays the species abundance data

if ~exist('hFig'),
    % Create figure window
    persistent hFig hLine
    hFig = figure(sim.figureNum);
%    set(hFig,'Name',['SimCom ' date()],'Position',[80,400,1100,450]);
    set(hFig,'Name',['SimCom ' date()],'Position',[474,630,802,314]);

    if sim.barPlot,
        bar([1:com.nSites],full(X'),'stacked');
        y_max = 10^log10(ceil(max(sum(X,1))));
        axis([0 com.nSites+1 0 y_max]);
        hLine = findobj('Type','patch');
    else
        plot([1:com.nSites],full(X'),'LineWidth',2);
        hold on
        plot([1, com.nSites],[1 1]*com.Xmin,'k--','LineWidth',2)
        hold off
        y_max = 10^full(log10(ceil(max(nonzeros(X)))));
        hLine = findobj('Type','Line');
        axis([0 com.nSites+1 0.9*com.Xmin 5]);
    end
    set(gca,'FontSize',20)
    xlabel('Spatial Position')
    ylabel('Species abundance')
    title(['t = ' num2str(t)])
    
else
    % Update figure
    if sim.barPlot,
        for s=1:com.nSpecies,
            % Updata the data for each of the lines
            data = full(X(s,:));
            yDat = [zeros(size(data)); data; data; zeros(size(data))];
            set(hLine(s),'YData',yDat);
        end
        y_max = 10^log10(ceil(max(sum(X,1))));
        axis([0 com.nSites+1 0 y_max]);
    else
        for s=1:com.nSpecies,
            % Updata the data for each of the lines
            set(hLine(s),'YData',full(X(s,:)))
        end
        set(get(gca,'Title'),'String',['t = ' num2str(t)])
    end
end

%drawnow
return
