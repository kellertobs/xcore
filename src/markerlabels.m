function markerlabels(y,lbl,FS)

ax = gca;

%-- Normalize y-data to [0,1]
if strcmp(ax.YScale,'log')
    yVals   = log10(y);
    yLimVals= log10(ax.YLim);
else
    yVals   = y;
    yLimVals= ax.YLim;
end
yNorm = (yVals - yLimVals(1)) / diff(yLimVals);

%-- Parameters
xNorm      = 1.02;      % just right of axis
minSpacing = 0.13;     % minimum vertical gap in normalized units

%-- Sort and cluster
[ys, order] = sort(yNorm);
clusters     = {};      % cell array of index vectors
current      = order(1);

for k = 2:numel(ys)
    if ys(k) - ys(k-1) < minSpacing
        current(end+1) = order(k);
    else
        clusters{end+1} = current; 
        current        = order(k);
    end
end
clusters{end+1} = current;

%-- Compute adjusted positions
yAdj = yNorm;  % initialize
for c = clusters
    idx = c{1};
    if numel(idx) > 1
        center = mean(yNorm(idx));
        n      = numel(idx);
        offsets = ((1:n) - (n+1)/2) * minSpacing;
        yAdj(idx) = center - offsets;
    end
    % singletons stay at yNorm(idx)
end

%-- Draw labels
for k = 1:numel(y)
    text(xNorm, yAdj(k), lbl{k}, FS{:}, ...
         'Units','normalized', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','middle', ...
         'Interpreter','latex');
end

end