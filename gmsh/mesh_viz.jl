using FerriteGmsh
using FerriteViz
import GLMakie

grid = togrid("grating.msh")
wireframe(grid, marksize=14, strokewidth=20, fontsize=25, nodelabel=true, celllabels=true)