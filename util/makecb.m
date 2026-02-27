% Create nice-looking colorbar
function cb = makecb(cminmax,cm,titlestr,sz,ax)


cb = colorbar;
clim(cminmax)
cb.TickLabelInterpreter = 'latex';
cb.FontSize = sz;
%drawnow;
cb.Label.String = titlestr;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = sz;
%ylabel(cb,titlestr,'Interpreter','latex','FontSize',sz);

if nargin == 4 
    colormap(cm)
elseif nargin > 4
    colormap(ax,cm)
end
