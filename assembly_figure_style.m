function S=assembly_figure_style()
%ASSEMBLY_FIGURE_STYLE Shared visual specification for assembly figures.
% Matches the cross-morph model-comparison figures.
S.morphOrder=["surface","molino","pachon"];
S.morphLabels=["Surface","Molino","Pachon"];
S.morphColors=[.85 .20 .20;.20 .65 .30;.95 .72 .10];
S.barWidth=.46;S.barAlpha=.78;S.barEdgeWidth=.8;
S.errorLineWidth=1.7;S.errorCapSize=8;
S.pointSize=25;S.pointAlpha=.82;S.pointEdgeWidth=.45;S.jitterWidth=.11;
S.curveLineWidth=2.5;S.fishLineWidth=.6;S.bandAlpha=.13;
S.referenceColor=[.65 .65 .65];S.referenceLineWidth=1;
S.pngResolution=300;
end
