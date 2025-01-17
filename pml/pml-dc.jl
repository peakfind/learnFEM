#=
# pml-dc.jl
#
# Use PML to truncate the unbounded domain.
# 
# Author: Jiayi Zhang
# Date:   17/01/2025
=#

using StaticArrays, Gmsh

"""
    setup_grid(d, n, b, lc)
    setup_grid(d, n, b, δ, lc, lp)

Generate the mesh by julia interface for Gmsh

# Arguments

- `d`: the period of the periodic grating
- `n`: the number of points on the profile of the grating
- `b`: the height of th truncated domain
- `lc`: target mesh size close to the point

- `δ`: the width of the pml layer
- `lp`: target mesh size in the pml layer

# Points
```
n+4 ------------- n+3
 |                 |
n+2 ------------- n+1
 |                 |
 1  -------------  n
```

# Lines
```
 .  ---- n+4 ----  .
n+5               n+3
 .  ---- n+1 ----  .
n+2                n 
 .  -- 1 ~ n-1 --  .
```
"""
function setup_grid(d, n, b, lc=0.5)
    h = d/(n-1)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    # Use a polyline to approximate the profile
    points = @MVector zeros(Int32, n+2)

    points[1] = gmsh.model.geo.addPoint(0, 0, 0, lc)
    points[n] = gmsh.model.geo.addPoint(d, 0, 0, lc)

    for i in 2:n-1
        points[i] = gmsh.model.geo.addPoint((i-1)*h, sin((i-1)*h), 0, lc) 
    end

    points[end-1] = gmsh.model.geo.addPoint(d, b, 0, lc) 
    points[end] = gmsh.model.geo.addPoint(0, b, 0, lc)

    # Add the lines
    lines = @MVector zeros(Int32, n+5)

    for i in 1:n+1
        lines[i] = gmsh.model.geo.addLine(points[i], points[i+1])
    end

    lines[end] = gmsh.model.geo.addLine(points[end], points[1])

    # Add the loop
    loop = gmsh.model.geo.addCurveLoop(lines)

    # Add the surface
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Generate the vtk file
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    # gmsh.write("pml.msh")
    gmsh.write("pml.vtk")
    gmsh.finalize()
end

function setup_grid(d, n, b, δ, lc=0.5, lp=0.1)
    h = d/(n-1)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    points = @MVector zeros(Int32, n+4)
    
    points[1] = gmsh.model.geo.addPoint(0, 0, 0, lc)
    points[n] = gmsh.model.geo.addPoint(d, 0, 0, lc)
    
    for i in 2:n-1
        points[i] = gmsh.model.geo.addPoint((i-1)*h, sin((i-1)*h), 0, lc) 
    end

    points[n+1] = gmsh.model.geo.addPoint(d, b, 0, lp)
    points[n+2] = gmsh.model.geo.addPoint(0, b, 0, lp)
    points[n+3] = gmsh.model.geo.addPoint(d, b+δ, 0, lp)
    points[n+4] = gmsh.model.geo.addPoint(0, b+δ, 0, lp)

    # Add the lines
    lines = @MVector zeros(Int32, n+5)

    for i in 1:n+1
        lines[i] = gmsh.model.geo.addLine(points[i], points[i+1])
    end

    lines[n+2] = gmsh.model.geo.addLine(points[n+2], points[1])
    lines[n+3] = gmsh.model.geo.addLine(points[n+1], points[n+3])
    lines[n+4] = gmsh.model.geo.addLine(points[n+3], points[n+4])
    lines[n+5] = gmsh.model.geo.addLine(points[n+4], points[n+2])

    # Add the loop
    lp1 = gmsh.model.geo.addCurveLoop(lines[1:n+2])
    lp2 = gmsh.model.geo.addCurveLoop([-lines[n+1], lines[n+3], lines[n+4], lines[n+5]])

    # Add the surface
    domain = gmsh.model.geo.addPlaneSurface([lp1])
    pml = gmsh.model.geo.addPlaneSurface([lp2])

    # Generate the vtk file
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    # gmsh.write("pml.msh")
    gmsh.write("pml.vtk")
    gmsh.finalize()
end

# d=2*pi, n=10, b=2
# setup_grid(2*pi, 10, 2)

# d=2*pi, n=10, b=2, δ=0.5
setup_grid(2*pi, 50, 2, 0.5)