function statsAxesClearAll(~,~,fig)
%statsAxesClearAll

% Reset all axes
f0 = get(fig.ax.scatter(1),'Children');
delete(f0); legend(fig.ax.scatter(1),'off');

f0 = get(fig.ax.conf(1),'Children');
delete(f0); legend(fig.ax.conf(1),'off');

f0 = get(fig.ax.load(1),'Children');
delete(f0); legend(fig.ax.load(1),'off');

f0 = get(fig.ax.spec(1),'Children');
delete(f0); legend(fig.ax.spec(1),'off');



end

