# Electromagnetic shower simulation

The objectives of the project were:

* Simulation of an electromagnetic shower initiated by a high energy photon produced in an aluminum material (graphite, iron and lead materials were also studied). 

* Compare the number of particles produced by photons of different energies.

* Compare the number of particles produced with the expected theoretical value.


Main parts of the code:

* Particle.cpp: class 'Particle' is defined. The most relevant characteristics of the particle are written here, such as energy, momentum and position.

* Formula.cpp: class 'Formula' is defined. The algorithms for computing the angles of the particles after being formed in 3 possible interactions are written here, as well as the calculations of the probabilities of their occurrence, and the energy lost due to inelastic collisions in the material.

* Propagator.cpp: class 'Propagator' is defined here and where it is made the propagation of each particle after the occurrence of each one of the 3 possible interactions: pair production, energy radiation by electron/positron and positron annihilation.

* 3ddrawing.cpp: draws the graphic representation of the electromagnetic shower in 3D. The trajectories of the particles that are formed are drawn until they reach a minimal energy of 50eV. From that point on they stop being considered in the simulation.

* histogram.cpp: draws the histograms for analysis. 

