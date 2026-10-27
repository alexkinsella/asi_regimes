function export_figure(fig, fig_dir, name)
%EXPORT_FIGURE Save a figure as a 300 dpi PNG in fig_dir.
%
% export_figure(fig, fig_dir, name) writes fig_dir/name.png, creating
% fig_dir if needed.
%
% Alex Kinsella, October 2026

if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
exportgraphics(fig, fullfile(fig_dir, [name '.png']), 'Resolution', 300)
end
