# Todo
## MUST be Done
- Investigate numberOfInterpolationSteps and paraviewTimestep crap - why/what...
- Scenario files should contain timestep size stuff
- Pass cellEdgeLength, Threadcount and C as parameters to reconstruction
- Reconstruction should save either reconstruction_implicit or reconstruction_explicit instead of just reconstruction
- Pass threadcount to simulation -> Do benchmarks to see how much thread creating/deletion costs

## SHOULD be Done
- Improve reconstruction performance (look for cache misses, memory allocation, unnecessary copying or work done, ...)
- Improve simulation performance (THREAD-POOLING!, look for cache misses, memory allocation, unnecessary copying or work done, ...)

## MAY be Done
- Implement quick and small particle renderer to view simulation result faster than in paraview and integrate automated video generation

## Thoughts for a Cleanup/Rewrite
- Single struct stat holds progress information (threads maybe, time remaining, saved frames, intermediate frames, async io, ...)
- Use Threadpool -> Performance + Threads that finished can help other threads
- Scenario "Inheritance" -> reuse default file for values that mostly never change (velocity factor, ETA etc.)
- Write down what data gets accessed when -> rearrange and optimize for that
- Redo timing system! -> Maybe talk about that in the next meeting
- Precompute adhesion and cohesion kernel if possible
- Check how low you can go with kernel precomputation granularity until it gets too inprecise -> MEASURE performance impact
- Consider variable amount of fluid particles (added at runtime -> emitter -> with different fluid mass -> HOW TO DO IT, next meeting maybe?)
- Consider moving boundary particles
- Pack all these data arrays in a struct or something to minimize parameter lists
- Method: Start from "scratch" -> MEASURE EVERY CHANGE FOR PERFORMANCE AND KEEP A HISTORY OR LOG OF IT!!!
