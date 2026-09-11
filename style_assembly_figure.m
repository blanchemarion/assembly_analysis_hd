function style_assembly_figure(fig)
%STYLE_ASSEMBLY_FIGURE Apply shared styling to every Cartesian/polar axes.
set(fig,'Color','w');
axesHandles=findall(fig,'Type','axes');
for k=1:numel(axesHandles)
 try,style_assembly_axes(axesHandles(k));catch,end
end
end
