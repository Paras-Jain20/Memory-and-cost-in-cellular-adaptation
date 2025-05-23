function cmap = color_map_func(varargin)
% Define the transition points
if(~isempty(varargin))
    % Define the number of colors in the colormap
    n = varargin{2};
    % Create a custom colormap with red, white, and blue
    cmap = zeros(n, 3);
    mid = varargin{1};
else
    % Define the number of colors in the colormap
    n = 256;
    % Create a custom colormap with red, white, and blue
    cmap = zeros(n, 3);
    mid = round(n/2);
end
% Red to white transition
for i = 1:mid
    cmap(i, 1) = 1; % Red channel
    cmap(i, 2) = (i - 1) / (mid - 1); % Green channel
    cmap(i, 3) = (i - 1) / (mid - 1); % Blue channel
end

% White to blue transition
for i = mid+1:n
    cmap(i, 1) = 1 - (i - mid) / (n - mid); % Red channel
    cmap(i, 2) = 1 - (i - mid) / (n - mid); % Green channel
    cmap(i, 3) = 1; % Blue channel remains constant
end
end