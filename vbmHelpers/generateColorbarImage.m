function generateColorbarImage(colormapName, numBins, limits, inv)
    % Generates and saves an image of the colorbar for a given colormap with specific bins
    % colormapName: any colormap name: see brainMontagerWithAtlas
    % numBins: number of bins in the colormap
    % limits: colorbar limits
    % inv: if true, colorbar will be inverted: default false

    cmap = colormap(colormapName);
    cmap = cmap(round(linspace(1, size(cmap, 1), numBins)), :);
    if inv
        cmap = flipud(cmap);
    end
    
    figure;
    colormap(cmap);c = colorbar;    
    c.Ticks = linspace(0, 1, numBins);
    numTicks = min(numBins, 10); % we need to limit this somehow, otherwise can be incomprehensible
    c.TickLabels = round(linspace(limits(1), limits(2), numTicks), 2); % Round to 2 decimal places for better readability
    c.Ticks = linspace(0, 1, numTicks);
    axis off; set(gcf, 'Color', 'white');
end
