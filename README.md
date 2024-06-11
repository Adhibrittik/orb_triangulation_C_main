# orb_triangulation_tetrahedral
Utilities concerning orbifold triangulations of finite covers of H^3/PGL(2,O_3). 
We can use TestForCovers.py to check if there are orientation and triangulation preserving covers from one finite cover of H^3/PGL(2,O_3) to another, and if these covers exist, TestForCovers.py also gives insight on which cusp goes where. 
SigToSeq.py extracts an orbifold destination sequence/orbifold triangulation (lifted from that of H^3/PGL(2,O_3)) from a Regina triangulation of a tetrahedral manifold consisting of regular ideal tetrahedra. Here, a tetrahedral manifold is a hyperbolic 3-manifold which has a decomposition into regular ideal tetrahedra (see the Fominykh-Garoufalidis-Goerner-Tarkaev-Vesnin paper "https://www.tandfonline.com/doi/abs/10.1080/10586458.2015.1114436"). 
SigToSeq.py needs to be run with regina-python (see "https://regina-normal.github.io/docs/man-regina-python.html"). 
