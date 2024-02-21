# Questions
- How to determine boundary particle rest density -> Use fluid rest density!
- How to add new fluid particles during the simulation => Different mass -> Just add particles with same mass
- Are there ways to couple parameters together (e.g. tension forces are mass dependent) -> NO!
- How much should rest density be increased to be stable -> Test yourself
- Report template? -> NO!

# Notes
- Report: 2-3 Pages -> SPH -> WCSPH -> PBF Explanation
- Report: 1-2 Pages -> Implementation

# Todo
- Check if the file(res to output res) directory copy in the CMakeLists file is necessary -> If not, remove res directory
- Update kernel_tests to test precomputated kernel against correct computation + general cleanup and revision
- Precompute adhesion and cohesion kernel if possible
- Precompute grad cubic function
- Use vtk input to reconstruction (saves time by not saving an extra file)
- Implement better method than manually resizing velocities, boundaryVolumes etc. after initial sampling
- Implement drag forces for fast particles (ref: exercise 2 - assignment 3c)
- Pass cellEdgeLength, Threadcount and C as parameters to reconstruction
- Reconstruction should save either reconstruction_implicit or reconstruction_explicit instead of just reconstruction
- Implement quick and small particle renderer to view simulation result faster than in paraview and integrate automated video generation
- Scenario "Inheritance" -> reuse default file for values that mostly never change (velocity factor, ETA etc.)
- Check how low you can go with kernel precomputation granularity until it gets too inprecise -> MEASURE performance impact
- Consider variable amount of fluid particles (added at runtime -> emitter -> with different fluid mass -> HOW TO DO IT, next meeting maybe?)
- Consider moving boundary particles
