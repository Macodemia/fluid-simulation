Task:
    Fix Parameters:
        Total Seconds: 5 seconds
        Eta: Used for computing h (multiplied with representative fluid particle diameter) -> reasonable: 1.2 -> h is smoothing length 
        Epsilon: Scaling parameter for velocity smoothing -> reasonable: 0.5
    Try Parameters:
        Particle Count (1000, 10000, 100000, 1000000)
        Max Time Step (0.001, 0.003, 0.008, 0.02, 0.08)
        Velocity Based Max Time Step Factor (0.25, 0.5, 0.75, 1, 2, 5)
        Number of Constraint Iterations -> (1, 2, 4, 8, 16) -> ONLY IMPLICIT: More Iterations: Closer rest density approximation. Converges very slowly, but low iteration counts produce visually good results
        Stiffness B: (1000, 2500, 5000, 7500, 10000) -> ONLY EXPLICIT: High Values: Instabilities, overshooting; Low Values: Particles go through boundary
        Height of fluid box -> (1, 1.5, 2, 4, 8)
    Comment on:
        Stability
        Compressibility
        Run Time

Test Scenes:
    Dam Break: Fluid box falls into box
    Double Dam Break: Two fluid boxes fall into box

Plan:
    1. Determine relationship between and reasonable values for MAX_TIMESTEP_SIZE and VELOCITY_TIME_STEP_FACTOR
        - Choose remaining parameters as safe/stable/reasonable as possible
        - MAX_TIMESTEP_SIZE = [0.002, 0.005, 0.008, 0.01, 0.02];
        - VELOCITY_TIME_STEP_FACTOR = [0.33, 0.66, 1, 2];
        - Particles: 10000
        - Fluidbox Height: 1
        - STIFFNESS_B: 1000 (As stable as possible, regarding explosions)
        - CONSTRAINT_ITERATIONS: 10 (As stable as possible)
        - Particle Radius: 0.0232079
        - Results:
              EXPLICIT:
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight1_particles10000: Stable and looking good 
            - gen_explicit_maxT0.002_velT0.66_stiff1000_fheight1_particles10000: Stable and looking good
            - gen_explicit_maxT0.002_velT1_stiff1000_fheight1_particles10000: Stable and looking good
            - gen_explicit_maxT0.002_velT2_stiff1000_fheight1_particles10000: Stable and looking good
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight1_particles10000: Stable but some pulsing
            - gen_explicit_maxT0.005_velT0.66_stiff1000_fheight1_particles10000: Extreme pulsing, some particles falling outside the box
            - gen_explicit_maxT0.005_velT1_stiff1000_fheight1_particles10000: Extreme pulsing, some particle falling outside the box
            - gen_explicit_maxT0.005_velT2_stiff1000_fheight1_particles10000: Constant extreme pulsing, many particles falling outside the box
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight1_particles10000: Stable but some pulsing
            - gen_explicit_maxT0.008_velT0.66_stiff1000_fheight1_particles10000: max time step too large; skip 
            - gen_explicit_maxT0.008_velT1_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.008_velT2_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.01_velT0.33_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.01_velT0.66_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.01_velT1_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.01_velT2_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.02_velT0.33_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.02_velT0.66_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.02_velT1_stiff1000_fheight1_particles10000: max time step too large; skip
            - gen_explicit_maxT0.02_velT2_stiff1000_fheight1_particles10000: max time step too large; skip
              IMPLICIT
            - gen_implicit_maxT0.002_velT0.33_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.002_velT0.66_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.002_velT1_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.002_velT2_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.005_velT0.33_iters10_fheight1_particles10000: Has one big bump/jump at the wall
            - gen_implicit_maxT0.005_velT0.66_iters10_fheight1_particles10000: Has two big bumps/jumps at first and second wall hit
            - gen_implicit_maxT0.005_velT1_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.005_velT2_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.008_velT0.33_iters10_fheight1_particles10000: Has one big and one small bump/jump
            - gen_implicit_maxT0.008_velT0.66_iters10_fheight1_particles10000: Has one big bump/jump
            - gen_implicit_maxT0.008_velT1_iters10_fheight1_particles10000: Has one big bump/jump
            - gen_implicit_maxT0.008_velT2_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.01_velT0.33_iters10_fheight1_particles10000: Has two small bumps/jumps
            - gen_implicit_maxT0.01_velT0.66_iters10_fheight1_particles10000: Has one very small bump/jump
            - gen_implicit_maxT0.01_velT1_iters10_fheight1_particles10000: Jumping/Pulsating at first ground impact
            - gen_implicit_maxT0.01_velT2_iters10_fheight1_particles10000: Stable, looking good but looses energy pretty fast
            - gen_implicit_maxT0.02_velT0.33_iters10_fheight1_particles10000: Multiple small/medium bumps/jumps
            - gen_implicit_maxT0.02_velT0.66_iters10_fheight1_particles10000: One small and one big bump/jump
            - gen_implicit_maxT0.02_velT1_iters10_fheight1_particles10000: Two big and some small bumps/jumps
            - gen_implicit_maxT0.02_velT2_iters10_fheight1_particles10000: Big jumps/bumps and particles falling outside the box
        - Findings:
            - Energy loss does not seem to depend on timestep size -> Test further and confirm/contradict
            - Max Timestep 0.02 (20ms) should not be used for either
            - Max Timestep 0.002 (2ms) both look good regardless of velocity factor as the max timestep is already small enough
            - Max Timestep should be below 0.005 (5ms) for explicit
            - What exact max time step and velocity timestep factor is optimal is still unclear. Here the max 0.002 timestep scenarios were always capped by this and not by the velocity factor. -> INVESTIGATE
            - Implicit may be able to use 0.01 (10ms) as timestep with high velocity timestep factors. In the 0.01 and velT2 case, the timestep was always capped to 0.01.
            - Implicit seems to not like variable time steps???? -> INVESTIGATE!
    2. Find out if implicit really does not like variable time steps
        - Disable velocity time step factor and retry previously (in step 1) unstable scenarios
        - Results:
            - gen_implicit_maxT0.005_velTxxx_iters10_fheight1_particles10000: looks good 
            - gen_implicit_maxT0.008_velTxxx_iters10_fheight1_particles10000: looks good 
            - gen_implicit_maxT0.01_velTxxx_iters10_fheight1_particles10000: looks good 
            - gen_implicit_maxT0.02_velTxxx_iters10_fheight1_particles10000: Still bad
        - Findings:
            - implicit does not like variable timesteps somehow, why?? Maybe that is a bug in our implementation??
    3. Test stiffness/constraint iterations for different fluid heights and relatively safe timesteps
        - Particles: 10000
        - Velocity Timestep Factor: 0.33 (COMPLETELY DISABLED FOR IMPLICIT)
        - Max Timesteps: 0.002, 0.005, 0.008
        - Fluid heights: 1, 1.5, 2
        - Stiffness B (explicit): 1000, 2500, 5000, 10000
        - Constraint Iterations (implicit): 2, 4, 8, 12
        - Runtimes ("not looked at" is replaced by similar result):
            EXPLCIIT:
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight1_particles10000:	296.892s	Looks good, little bouncing
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight1_particles10000:	299.313s	Looks ok, some bouncing
            - gen_explicit_maxT0.008_velT0.33_stiff2500_fheight1_particles10000:	337.826s	bad
            - gen_explicit_maxT0.005_velT0.33_stiff2500_fheight1_particles10000:	342.646s	Weird drop at the beginning, little bouncing throughout
            - gen_explicit_maxT0.008_velT0.33_stiff5000_fheight1_particles10000:	355.863s	bad
            - gen_explicit_maxT0.008_velT0.33_stiff10000_fheight1_particles10000:	358.244s	bad
            - gen_explicit_maxT0.005_velT0.33_stiff2500_fheight1.5_particles10000:	362.626s	Weird drop at the beginning, little bouncing throughout
            - gen_explicit_maxT0.005_velT0.33_stiff5000_fheight1_particles10000:	363.020s	Box explodes, big bouncing throughout
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight1.5_particles10000:	363.145s	Looks ok, first fall is a bit weird
            - gen_explicit_maxT0.005_velT0.33_stiff10000_fheight1_particles10000:	364.905s	Box explosion, constant bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight1.5_particles10000:	379.510s	Looks good, one particle falling through
            - gen_explicit_maxT0.005_velT0.33_stiff5000_fheight1.5_particles10000:	397.717s	Box explodes, big bouncing throughout
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight2_particles10000:	397.089s	Looks ok, tiny bounces throughout
            - gen_explicit_maxT0.005_velT0.33_stiff10000_fheight1.5_particles10000:	402.591s	Box explosion, constant bouncing
            - gen_explicit_maxT0.008_velT0.33_stiff2500_fheight1.5_particles10000:	404.418s	bad
            - gen_explicit_maxT0.008_velT0.33_stiff10000_fheight1.5_particles10000:	408.169s	bad
            - gen_explicit_maxT0.008_velT0.33_stiff5000_fheight1.5_particles10000:	410.166s	bad
            - gen_explicit_maxT0.005_velT0.33_stiff10000_fheight2_particles10000:	414.892s	Box explosion, constant bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight2_particles10000:	425.662s	Looks good, little bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff5000_fheight2_particles10000:	428.074s	Box explodes, big bouncing throughout
            - gen_explicit_maxT0.008_velT0.33_stiff2500_fheight2_particles10000:	429.727s	bad
            - gen_explicit_maxT0.002_velT0.33_stiff2500_fheight1_particles10000:	430.002s	looked good
            - gen_explicit_maxT0.005_velT0.33_stiff2500_fheight2_particles10000:	430.881s	Particles falling through, medium bouncing throughout 
            - gen_explicit_maxT0.008_velT0.33_stiff10000_fheight2_particles10000:	434.332s	bad
            - gen_explicit_maxT0.008_velT0.33_stiff5000_fheight2_particles10000:	434.755s	bad
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight1_particles10000:	439.796s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff5000_fheight1_particles10000:	454.648s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff10000_fheight1_particles10000:	463.333s	Box explodes, multiple bounces
            - gen_explicit_maxT0.002_velT0.33_stiff5000_fheight1.5_particles10000:	521.333s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight1.5_particles10000:	522.063s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff10000_fheight1.5_particles10000:	528.544s	Box slightly explodes, Multiple bounces
            - gen_explicit_maxT0.002_velT0.33_stiff2500_fheight1.5_particles10000:	537.043s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff5000_fheight2_particles10000:	614.667s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight2_particles10000:	631.015s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff2500_fheight2_particles10000:	632.311s	looked good
            - gen_explicit_maxT0.002_velT0.33_stiff10000_fheight2_particles10000:	647.900s	Box slightly explodes, multiple bounces
            IMPLICIT:
            - gen_implicit_maxT0.008_velT0.33_iters2_fheight1_particles10000:		245.726s	some particles fall trough, tiny tiny jitter but looks good
            - gen_implicit_maxT0.008_velT0.33_iters4_fheight1_particles10000:		251.509s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters8_fheight1_particles10000:		273.960s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters2_fheight1.5_particles10000:		289.241s	looks good, tiny tiny jitter
            - gen_implicit_maxT0.008_velT0.33_iters4_fheight1.5_particles10000:		302.883s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters12_fheight1_particles10000:		303.215s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters8_fheight1.5_particles10000:		340.420s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters12_fheight1.5_particles10000:	349.423s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters2_fheight2_particles10000:		351.796s	worst case for maxT0.008, looks ok but some jitter -> worse than best case
            - gen_implicit_maxT0.005_velT0.33_iters2_fheight1_particles10000:		359.936s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters4_fheight2_particles10000:		363.486s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters4_fheight1_particles10000:		379.528s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters8_fheight2_particles10000:		391.448s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters8_fheight1_particles10000:		409.666s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters2_fheight1.5_particles10000:		419.539s	looks good
            - gen_implicit_maxT0.008_velT0.33_iters12_fheight2_particles10000:		431.636s	best case for maxT0.008, looks good -> better than worst case
            - gen_implicit_maxT0.005_velT0.33_iters4_fheight1.5_particles10000:		441.320s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters12_fheight1_particles10000:		443.592s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters8_fheight1.5_particles10000:		497.801s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters2_fheight2_particles10000:		513.669s	looks good (worst case for maxT0.005 iters2), looks "the same" for best case maxT0.005 iters12
            - gen_implicit_maxT0.005_velT0.33_iters12_fheight1.5_particles10000:	521.926s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters4_fheight2_particles10000:		538.836s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters8_fheight2_particles10000:		583.983s	looks good
            - gen_implicit_maxT0.005_velT0.33_iters12_fheight2_particles10000:		606.219s	looks good (best case maxT0.005 iters12), looks "the same" for worst case for maxT0.005 iters2
            - gen_implicit_maxT0.002_velT0.33_iters2_fheight1_particles10000:		877.072s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters4_fheight1_particles10000:		911.285s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters8_fheight1_particles10000:		1021.12s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters2_fheight1.5_particles10000:		1036.20s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters4_fheight1.5_particles10000:		1084.13s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters12_fheight1_particles10000:		1096.52s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters8_fheight1.5_particles10000:		1198.65s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters2_fheight2_particles10000:		1257.25s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters12_fheight1.5_particles10000:	1267.08s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters4_fheight2_particles10000:		1319.59s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters8_fheight2_particles10000:		1437.98s	looks good
            - gen_implicit_maxT0.002_velT0.33_iters12_fheight2_particles10000:		1537.54s	looks good
        - Results (unsorted and original filetree order):
              EXPLICIT:
            - gen_explicit_maxT0.002_velT0.33_stiff10000_fheight1.5_particles10000: Box slightly explodes, Multiple bounces
            - gen_explicit_maxT0.002_velT0.33_stiff10000_fheight1_particles10000: Box explodes, multiple bounces
            - gen_explicit_maxT0.002_velT0.33_stiff10000_fheight2_particles10000: Box slightly explodes, multiple bounces
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight1.5_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight1_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff1000_fheight2_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff2500_fheight1.5_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff2500_fheight1_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff2500_fheight2_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff5000_fheight1.5_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff5000_fheight1_particles10000: Looks good
            - gen_explicit_maxT0.002_velT0.33_stiff5000_fheight2_particles10000: Looks good
            - gen_explicit_maxT0.005_velT0.33_stiff10000_fheight1.5_particles10000: Box explosion, constant bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff10000_fheight1_particles10000: Box explosion, constant bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff10000_fheight2_particles10000: Box explosion, constant bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight1.5_particles10000: Looks good, one particle falling through
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight1_particles10000: Looks good, little bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff1000_fheight2_particles10000: Looks good, little bouncing
            - gen_explicit_maxT0.005_velT0.33_stiff2500_fheight1.5_particles10000: Weird drop at the beginning, little bouncing throughout
            - gen_explicit_maxT0.005_velT0.33_stiff2500_fheight1_particles10000: Weird drop at the beginning, little bouncing throughout
            - gen_explicit_maxT0.005_velT0.33_stiff2500_fheight2_particles10000: Particles falling through, medium bouncing throughout
            - gen_explicit_maxT0.005_velT0.33_stiff5000_fheight1.5_particles10000: Box explodes, big bouncing throughout
            - gen_explicit_maxT0.005_velT0.33_stiff5000_fheight1_particles10000: Box explores, big bouncing throughout 
            - gen_explicit_maxT0.005_velT0.33_stiff5000_fheight2_particles10000: Box explores, big bouncing throughout
            - gen_explicit_maxT0.008_velT0.33_stiff10000_fheight1.5_particles10000: not looked at, "easier" scenarios are already bad 
            - gen_explicit_maxT0.008_velT0.33_stiff10000_fheight1_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff10000_fheight2_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight1.5_particles10000: Looks ok, first fall is a bit weird
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight1_particles10000: Looks ok, some bouncing
            - gen_explicit_maxT0.008_velT0.33_stiff1000_fheight2_particles10000: Looks ok, tiny bounces throughout
            - gen_explicit_maxT0.008_velT0.33_stiff2500_fheight1.5_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff2500_fheight1_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff2500_fheight2_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff5000_fheight1.5_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff5000_fheight1_particles10000: not looked at, "easier" scenarios are already bad
            - gen_explicit_maxT0.008_velT0.33_stiff5000_fheight2_particles10000: not looked at, "easier" scenarios are already bad
              IMPLICIT:
            - gen_implicit_maxT0.002_velT0.33_iters12_fheight1.5_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters12_fheight1_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters12_fheight2_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters2_fheight1.5_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters2_fheight1_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters2_fheight2_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters4_fheight1.5_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters4_fheight1_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters4_fheight2_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters8_fheight1.5_particles10000: looks good
            - gen_implicit_maxT0.002_velT0.33_iters8_fheight1_particles10000: not looked at, as maxT0.002 is fine for implicit
            - gen_implicit_maxT0.002_velT0.33_iters8_fheight2_particles10000: not looked at, as maxT0.002 is fine for implicit
            - gen_implicit_maxT0.005_velT0.33_iters12_fheight1.5_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters12_fheight1_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters12_fheight2_particles10000: looks good (best case maxT0.005 iters12), looks "the same" for worst case for maxT0.005 iters2
            - gen_implicit_maxT0.005_velT0.33_iters2_fheight1.5_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters2_fheight1_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters2_fheight2_particles10000: looks good (worst case for maxT0.005 iters2), looks "the same" for best case maxT0.005 iters12
            - gen_implicit_maxT0.005_velT0.33_iters4_fheight1.5_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters4_fheight1_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters4_fheight2_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters8_fheight1.5_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters8_fheight1_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.005_velT0.33_iters8_fheight2_particles10000: not looked at as worst and best case for maxT0.005 are almost identical
            - gen_implicit_maxT0.008_velT0.33_iters12_fheight1.5_particles10000: not looked at, as fheight2 is already looking good
            - gen_implicit_maxT0.008_velT0.33_iters12_fheight1_particles10000: not looked at, as fheight2 is already looking good
            - gen_implicit_maxT0.008_velT0.33_iters12_fheight2_particles10000: best case for maxT0.008, looks good -> better than worst case
            - gen_implicit_maxT0.008_velT0.33_iters2_fheight1.5_particles10000: looks good, tiny tiny jitter
            - gen_implicit_maxT0.008_velT0.33_iters2_fheight1_particles10000: some particles fall trough, tiny tiny jitter but looks good
            - gen_implicit_maxT0.008_velT0.33_iters2_fheight2_particles10000: worst case for maxT0.008, looks ok but some jitter -> worse than best case
            - gen_implicit_maxT0.008_velT0.33_iters4_fheight1.5_particles10000: looks good
            - gen_implicit_maxT0.008_velT0.33_iters4_fheight1_particles10000: looks good
            - gen_implicit_maxT0.008_velT0.33_iters4_fheight2_particles10000: looks good
            - gen_implicit_maxT0.008_velT0.33_iters8_fheight1.5_particles10000: looks good
            - gen_implicit_maxT0.008_velT0.33_iters8_fheight1_particles10000: looks good
            - gen_implicit_maxT0.008_velT0.33_iters8_fheight2_particles10000: looks good
        - Findings:
            - Below 0.008 (8ms), iteration count does not really matter for implicit
            - At 0.008 at least 4 iters should be used
            - Explicit cannot handle stiffness of 10000
            - Explicit cannot handle stiffness of 5000 above or equal to maxT0.005
            - Seems reasonable: maxT0.002 and stiffness smaller or equal 5000
            - Larger timesteps seem to loose very slightly less energy (observed at 0.002 and 0.008) 
    4. Try to get "the same" result for explicit and implicit (dam break scene)
        - Explicit Parameters:
        - Implicit Parameters:

TODO: Find optimum for iteration count / max time step (result - run time wise)
TODO: We only know that, for explicit, 0.002 is fine and 0.005 is too large -> get more fine grained analysis
TODO: Maybe remove 0.002 max timestep from further tests as it is stable for both (maybe only use for run time calculation or something)
TODO: Try different number of particles
