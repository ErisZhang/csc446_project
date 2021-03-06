
### todo

1. c++ impl 
    + minimize allocation, i.e. B matrix
    + solver
        + `ldlt` -> `llt` (stable on positive definite matrix) 
        + set iteration number to less, monitor residue ...


## iterative methods 


## Gradient estimation on discrete mesh 

+ [2018_gradient_field_estimation_on_triangle_meshes](2018_gradient_field_estimation_on_triangle_meshes.pdf)
    + estimate the gradient operator for discrete mesh,
    + reference for 
        + https://github.com/libigl/libigl/blob/508cb9940f4d1e8e54137d5afe2fd2eb9c4dc672/include/igl/grad.h
        + https://github.com/alecjacobson/gptoolbox/blob/96783f4e17a3e8c26a9480c31dbcdafab5629afb/mesh/grad.m


## Wiki

+ constitutive equation
    + https://en.wikipedia.org/wiki/Constitutive_equation
    + https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Constitutive_relations#Isotropic_materials

+ stress
    + https://en.wikipedia.org/wiki/Stress_(mechanics)

+ linear elasticity
    + https://en.wikipedia.org/wiki/Linear_elasticity

+ direct stiffness method
    + https://en.wikipedia.org/wiki/Direct_stiffness_method


+ MIT 16.20 structural mechanics
    + http://web.mit.edu/16.20/homepage/index.html
    + super good introduction to continuum mechanics 
    + some points that stands out
        + good examples
        + usage of tensor notation 


## textbooks / chapters / slides / reports 

+ [6.2_iterative_methods](6.2_iterative_methods.pdf)
    + convergence tst
        + M = I - P^{-1}A must have \rho(M) < 1
        + 0.9 0.99 is ok

+ [text_2003_iterative_methods_for_sparse_linear_systems2ed](text_2003_iterative_methods_for_sparse_linear_systems2ed.pdf)
    + sparse linear solvers, preconditioning, convergence proofs ...


+ [text_2009_structural_analysis_with_applications_to_aerospace_structures](text_2009_structural_analysis_with_applications_to_aerospace_structures.pdf)
    + great book on introducing basic concepts in linear elasticity theory
    + chapter 1-4, skip special case analytic solution procedures
    + but did not provide enough detail on FEM to solve the stress field

+ [text_2002_the_linearized_theory_of_elasticity](text_2002_the_linearized_theory_of_elasticity.pdf)
    + pretty math intensive, but quite expressive and clean (tensor notation mostly)
    + didn't really talk about how to solve linear elasticity equation using FEM

+ [text_2010_the_finite_element_method_theory_implementation_and_practice](text_2010_the_finite_element_method_theory_implementation_and_practice.pdf)
    + details on FEM
    + chapter 11 specifically on FEM on linear elasticity problems
    + accompanying assignment: http://www.it.uu.se/edu/course/homepage/finmet2/vt14/material/proj1.pdf

+ siggraph 2012 course in FEM and structural analysis [part1](siggraph_2012_part1_classical_fem_method_and_discretization_methodology.pdf) and [part2](siggraph_2012_part2_model_reduction.pdf)

+ online tutorial http://www.continuummechanics.org/index.html
    + very good !
    + intro to tensor notations


+ [course_linear_elasticity_variational_formulation](course_linear_elasticity_variational_formulation.pdf)
    + simple formulation 

+ [ch4_numerical_integration_methods_to_evaluate_triple_integrals_using_generalized_gaussian_quadrature](ch4_numerical_integration_methods_to_evaluate_triple_integrals_using_generalized_gaussian_quadrature.pdf)
    + gaussian quadrature in 3d

+ [ch16_application_of_assembly_of_finite_element_methods_on_graphics_processors_for_real_time_elastodynamics]
    + good explanation of FEM method 
    + GPU algos/impls for FEM on GPU

+ [08_finite_elements_basisfunctions](08_finite_elements_basisfunctions.pdf)
    + basis function for 1d/2d

+ [me582_ch3_fea_on_thermofluids](me582_ch3_fea_on_thermofluids.pdf)


+ how to apply dirichlet boundary constraints to global stiffness matrix
    + http://podgorskiy.com/spblog/304/writing-a-fem-solver-in-less-the-180-lines-of-code
    + http://solidmechanics.org/text/Chapter7_2/Chapter7_2.htm
    


## Softwares


+ Eigen seem to support matrix-free iterative solver with conjugate gradient 
    + https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html



## computer implementation for fem 

+ [libigl tutorial on laplace equation](https://libigl.github.io/tutorial/#laplace-equation)
    + min quad with fixed has details on usage of lagrange multiplier in case of linear constraints on `z`

+ [15_linear_tetrahedron](15_linear_tetrahedron.pdf)
    + formulation of fem on tetrahedron with linear basis function

+ online tutorial 
    + on tetrahedrons [part1](http://what-when-how.com/the-finite-element-method/fem-for-3d-solids-finite-element-method-part-1/)
        + derivation of barycentric coordinate, matrix form of stiffness operator, interpolation ...
    + on hexahedrons [part2](http://what-when-how.com/the-finite-element-method/fem-for-3d-solids-finite-element-method-part-2/)

+ [computational fabrication simulation slides](slides_compfab_05_simulation_ii.pdf)
    + pretty helpful on going from energy -> nodal displacement formulation for element-wise stiffness `K_e`

+ [2002_matlab_implementation_of_the_finite_element_method_in_elasticity](2002_matlab_implementation_of_the_finite_element_method_in_elasticity.pdf)
    + 2d linear elasticity in matlab
    + non-vectorized, but simple to understand
    + might be helpful to use the terminologies for writing the report
    + has 3d impl for tetrahedral elements

+ [2007_vectorized_matlab_codes_for_linear_2d_elasticity](2007_vectorized_matlab_codes_for_linear_2d_elasticity.pdf)
    + 2d linear elasticity in matlab 
        + gptoolbox's `linear_elasticity` based off of this


+ [2010_cuda_for_realtime_multigrid_finite_element_simulation_of_soft_tissue_deformations](2010_cuda_for_realtime_multigrid_finite_element_simulation_of_soft_tissue_deformations.pdf)
    + cuda implementation 
    + realtime performance for linear elasticity

+ [2011_a_realtime_multigrid_finite_hexahedra_method_for_elasticity_simulation_using_CUDA](2011_a_realtime_multigrid_finite_hexahedra_method_for_elasticity_simulation_using_CUDA.pdf)
    + cuda impl for hexa-hedral element
    + matrix-free


## FEM methods for linear elasticity

+ [2006_comparative_analysis_of_PCG_solvers_for_voxel_FEM_systems](2006_comparative_analysis_of_PCG_solvers_for_voxel_FEM_systems.pdf)
    + preconditioning methods for a conjugate gradient solver for an elliptic BVP 
    + example given is a generalization of poisson equation (not really related to linear elasticity)
    + but has some details on basis function in 3d

+ [2008_finite_element_methods_for_linear_elasticity](2008_finite_element_methods_for_linear_elasticity.pdf)
    + a mathematically rigorous formulation of several weak form of linear elasticity
    + concise!
        + various variational formulation

+ [2011_finite_element_analysis_for_linear_elastic_solids_based_on_subdivision_schemes](2011_finite_element_analysis_for_linear_elastic_solids_based_on_subdivision_schemes.pdf)
    + details on assembly stiffness `k` matrix, similar to CompFab A4
    + voxel based formulation of linear elasticity

+ [2014_large_scale_finite_element_analysis_via_assembly_free_deflated_conjugate_gradient](2014_large_scale_finite_element_analysis_via_assembly_free_deflated_conjugate_gradient.pdf)
    + idea
        + assembly free deflated conjugate gradient
    + a good survey of existing methods
        + preconditioning
        + multigrid/deflation
    + talked about how element congruency can accelerate assembly free system by 
        + reducing memory usage of storing/retrieving stiffness/deflation matrix
        + conjugate gradient an be accelerated in assembly free mode
    + element congruency
        + same element stiffness matrix!
    + results
        + 3.15M dof / 83,000 elements
            + Jacobi-PCG: 1741 iterations and 245seconds
            + element congruency: CPU 52s GPU 39s
    

+ [2016_multigrid_preconditioning_of_linear_elasticity_equations](2016_multigrid_preconditioning_of_linear_elasticity_equations.pdf)
    + good background on 
        + deriving variational formulation of linear elasticity
        + link to other references
    + more details on
        + iterative methods
        + multigrid methods
    + good references to 
        + 3d basis for virtual element method (https://arxiv.org/pdf/1311.0932.pdf)
        + 

+ [matrix_free_voxel_based_finite_element_method_for_materials_with_heterogeneous_microstructures](matrix_free_voxel_based_finite_element_method_for_materials_with_heterogeneous_microstructures.pdf)
    + a thesis that is quite well written
    + matrix-free algorithm on voxels ! exactly what we wanted ... base our project on this ...


## Multigrids


+ [a_multigrid_tutorial_2ed](a_multigrid_tutorial_2ed.pdf)


## graphics paper based on linear elasticity

+ [2012_stress_relief_improving_structural_strength_of_3D_printable_objects](2012_stress_relief_improving_structural_strength_of_3D_printable_objects.pdf)
    + abstract 
        + automatic detection and correction of problematic cases for 3D printing objects to make more "structurally sound" while retaining visual similarity
        + idea
            + find high structural stress positions
            + remediate with 
                + hollowing (weights)
                + thickening (increases strength of thin)
                + struct insertion (prevent nonrigid deformation)
        + steps
            + 
    + 1 related work/paper
        + paper on structurally stable method by modifying shape to minimize a certain objective function
    + 3 printability analysis
        + structural loads
            + permanent load (gravity)
                + given _default upright orientation_
                + and find additional orientations with some heuristic
            + imposed load (gravity + pinch grip)
                + 2 finger pinch (location determined with some heuristic)
        + structural analysis
            + problem: linear elasticity problem with a homogenous material
            + tet mesh 
            + discrete linear elasticity problem using FEM with quadratic tetrahedral elements
                + `Kd=F`
                    + `K` stiffness matrix
                    + `d` deformation of mesh 
                    + `F` force
                + use deformation to compute 
                    + strain: Cauchy linear strain tensor
                    + stress: hook's law
                + von mises yield criterion to determine structural problems
     


+ [2013_worst_case_structural_analysis](2013_worst_case_structural_analysis.pdf)
    + abstract
        + identify structural problems in 3D printing objects based on _geometry_ and _material property_ only, without assumptions on loads.
        + idea
            + constrained optimization to determine _worst case load distribution_ for a shape that will cause high local stress or large deformations (i.e. minimal force to break it or severely deform it)
            + approximate method to make it run faster
    + 2 related work / paper
        + 2012 paper finds structural weakness considering gravity loads and gripping shape with 2 fingers
        + 3d printing processes
            + materials
                + brittle: need to know where material likely to break
                + ductile: need to know where material likely to be plastic
            + goal:
                + find worst case loads (that leads to maximal damage)
                    + norm of stress
                    + maximal displacements
    + 3 worst case structural analysis
        + 



+ [2016_stochastic_structural_analysis_for_context_aware_design_and_fabrication](2016_stochastic_structural_analysis_for_context_aware_design_and_fabrication.pdf)
    + goal 
        + probability of failure as a measure of reliability
    + background
        + FEM
            + tool for failure analysis 
            + drawbacks
                + spatially varying ratings not intuitive
                + worst-case scenario may never happen (in real-world use)
        + stochastic FEM
            + idea
                + assume material property fixed, apply load stochastically with unknown variance
                    + i.e. interaction of linearly elastic, infinitesimally deforming objects with the world, both persistent and transient
                    + simulated with rigid body simulation
            + compute probability of failure
                + i.e. a toy can survive 99% of interactions in real-world
            + force distribution instead of single instance of applied force
            + method
                + method of generating realistic force distributions (rigid-body physics engine generates force samples on surface) (speed by 1 order of magnitude from MC based appr)
                + SVD -> low-dim reduced force-space
                + compute failure probability density function (probability that max stress experienced by an object will exceed a given threshold)
        + inverse design problem
            + given
                + mesh
                + volumetric material assignment
                + usage scenario
                + minimum fracture probability
            + automatically correct design flaws
                + reduce weight
                + guarantee probablity of breaking less than specified value
    + 2 related work/paper
        + FEM for (worst-case) structural analysis
            + `2012_stress_relief_improving_structural_strength_of_3D_printable_objects`
            + `2013_worst_case_structural_analysis`
            + weakness map ... areas of failure
        + SFEM
            + different kinds (monte carlo/ perturbation)
        + topology optimization
            + topology of an object is optimized to meet certain requirement
    + 3 background on structural analysis
        + FEM for linearly elastic objects
            + compute stress induced on an objct by the applied force
            + steps
                + discretization with mesh of volumetric cell
                + static equilibrium `Ku=f`
                    + `u` is nodal displacement
                + stress `sigma=CBu`
                    + `sigma` is cauchy stres
        + detecting failure using yield stress
            + fracture
                + `von Mises stress > yield stress (material specific)`
        + stochastic FEM
            + `u` and `f` replaced by random variables
    + 4 fracture probabilities using SFEM
        + generate force distributions
            + assumptions
                + object deforms infinitesimally
                + object instantaneously returns to its undeformed state after a force is applied
            + goodness
                + assume geometry and inertial property fixed over time (hence can be simulated as a rigid body)

