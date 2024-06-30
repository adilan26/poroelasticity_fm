[xi, yi] = meshgrid( min(x):1: max(x),  min(y):0.01: max(y));
zi = griddata(x,y, (variable_plot), xi,yi);
contourf(xi,yi,(zi),100,'ShowText','off','edgecolor','none');
c = colorbar;

xlabel ('{\textbf{Horizontal distance}} {\textbf{(m)}}  ','Interpreter', 'Latex' );
ylabel ('{\textbf{Vertical distance}} {\textbf{(m)}}  ','Interpreter', 'Latex' );

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontWeight', 'bold');