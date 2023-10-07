# BiharmonicNavier

* Title: A $C^0$ finite element method for the biharmonic problem with Navier boundary conditions in a polygonal domain.

* Abstract: In this paper, we study the biharmonic equation with Navier boundary conditions in a polygonal domain. In particular, we propose a method that effectively decouples the fourth-order problem as a system of Poisson equations. Our method differs from the naive mixed method that leads to two Poisson problems but only applies to convex domains; our decomposition involves a third Poisson equation to confine the solution in the correct function space and therefore can be used in both convex and nonconvex domains. A C0 finite element algorithm is in turn proposed to solve the resulting system. In addition, we derive optimal error estimates for the numerical solution on both quasi-uniform meshes and graded meshes. Numerical test results are presented to justify the theoretical findings.

* This code recovers the example in [1] in an L-shaped domain.

run
main

* The current code only shows the case with an L-shaped domain, if you need to calculate more general angles, or more assistance is needed, please reach out to PeimengYin@gmail.com

* To cite this code or the original article, please cite [1] below.

* [1] Hengguang Li, Peimeng Yin and Zhimin Zhang. A $C^0$ finite element method for the biharmonic problem with Navier boundary conditions in a polygonal domain. IMA Journal of Numerical Analysis, 43:1779-1801, 2023.

* The code is developed based on iFEM by Long Chen.
