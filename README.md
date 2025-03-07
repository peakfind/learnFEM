# learnFEM
 some julia codes for learning fem

 - ch1, ch2, and ch3 codes are based on [M.G. Larson & F. Bengzon, The Finite Element Method: Theory, Implementation, and Applications](https://link.springer.com/book/10.1007/978-3-642-33287-6).
 - gmsh contains some .geo files and some julia scripts to use [Gmsh](http://www.gmsh.info/).
 - pbc contains .jl files which implement Dirichlet Periodic boundary condition. `pbc1.jl` use `Ferrite` to generate the grid. `pbc2.jl` uses `Gmsh` to generate the grid.
 - Be careful! In `pml/periodic_surface.jl` we need to deal the `BCValues()` in `ConstraintHandler()`.
