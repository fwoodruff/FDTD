# FDTD
Finite Difference Time Domain solver for modelling Terahertz Waveplates

Imagine a region of space containing arbitrary conducting and dielectric materials.
A microwave pulse enters the region. How does the system evolve?

I approximated continuous space as a lattice and gave each point a value for conductivity σ, permittivity ε, and electromagnetic field {**E**, **B**}. I rearranged Maxwell's Equations to see how the fields at each point would evolve in a time dt, and then integrated over time. A numerical method is used to limit reflections at the boundaries of the modelled region.

There is some simple multithreading and a matplotlib animation output.
