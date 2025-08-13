function markerlabels(y,lbl,FS,side)

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
if strcmp(side,'right')
    xNorm      = 1.02;     % just right of axis
else
    xNorm      = 0.98;     % just left of axis
end
minSpacing = 0.2;     % minimum vertical gap in normalized units

%-- Sort and cluster
[ys, order] = sort(yNorm,'descend');
clusters     = {};      % cell array of index vectors
current      = order(1);
dist         = abs(diff(ys));

for k = 2:numel(ys)
    if dist(k-1) < minSpacing
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
        offsets = ((1:n) - (n+1)/2) * (minSpacing-0.07);
        yAdj(idx) = center - offsets;
    end
    % singletons stay at yNorm(idx)
end

%-- Draw labels
if strcmp(side,'right')
    align = 'left';
else
    align = 'right';
end
for k = 1:numel(y)
    text(xNorm, yAdj(k), lbl{k}, FS{:}, ...
         'Units','normalized', ...
         'HorizontalAlignment',align, ...
         'VerticalAlignment','middle', ...
         'Interpreter','latex');
end

end