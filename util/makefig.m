% Make good figure window with optional width and height
function fig = makefig(width, height)

    % Set default size if arguments are not given
    if nargin < 1
        width = 1800;
    end
    if nargin < 2
        height = 1800;
    end

    % Position format: [left bottom width height]
    fig = figure('Color', [1 1 1], 'Position', [100 100 width height]);
end