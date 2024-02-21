#!/usr/bin/env node
const fs = require('fs');
process.chdir(__dirname);

const numberOfParticles_values = [1000, 10000, 50000];
const MAX_TIMESTEP_SIZE_values = [0.0]; // CURRENT HACK: This value is different for explicit and implicit and is therefore set in the corresponding loop
const VELOCITY_TIME_STEP_FACTOR_values = [0.5];
const fluid_height_values = [2];

const STIFFNESS_B_values = [3000]; // Only explicit
const CONSTRAINT_ITERATIONS_values = [4]; // Only implicit

let explicit_combinations = [];
for (let numberOfParticles of numberOfParticles_values) {
    for (let MAX_TIMESTEP_SIZE of MAX_TIMESTEP_SIZE_values) {
        MAX_TIMESTEP_SIZE = 0.002; // TODO: HACK!
        for (let VELOCITY_TIME_STEP_FACTOR of VELOCITY_TIME_STEP_FACTOR_values) {
            for (let STIFFNESS_B of STIFFNESS_B_values) {
                for (let fluid_height of fluid_height_values) {
                    explicit_combinations.push({
                        'numberOfParticles': numberOfParticles,
                        'MAX_TIMESTEP_SIZE': MAX_TIMESTEP_SIZE,
                        'VELOCITY_TIME_STEP_FACTOR': VELOCITY_TIME_STEP_FACTOR,
                        'STIFFNESS_B': STIFFNESS_B,
                        'fluid_height': fluid_height
                    });
                }
            }
        }
    }
}

let implicit_combinations = [];
for (let numberOfParticles of numberOfParticles_values) {
    for (let MAX_TIMESTEP_SIZE of MAX_TIMESTEP_SIZE_values) {
        MAX_TIMESTEP_SIZE = 0.008; // TODO: HACK!
        for (let VELOCITY_TIME_STEP_FACTOR of VELOCITY_TIME_STEP_FACTOR_values) {
            for (let CONSTRAINT_ITERATIONS of CONSTRAINT_ITERATIONS_values) {
                for (let fluid_height of fluid_height_values) {
                    implicit_combinations.push({
                        'numberOfParticles': numberOfParticles,
                        'MAX_TIMESTEP_SIZE': MAX_TIMESTEP_SIZE,
                        'VELOCITY_TIME_STEP_FACTOR': VELOCITY_TIME_STEP_FACTOR,
                        'CONSTRAINT_ITERATIONS': CONSTRAINT_ITERATIONS,
                        'fluid_height': fluid_height
                    });
                }
            }
        }
    }
}

// TODO: The output scenario is currently hardcoded as "Dam Break"

console.log('Combination Count (Explicit): ' + explicit_combinations.length);
console.log('Combination Count (Implicit): ' + implicit_combinations.length);

// Explicit Scenario
for (const combination of explicit_combinations) {
    let explicitScenario = `\
TOTAL_SECONDS_TO_SIMULATE: 4.5
MAX_TIMESTEP_SIZE: ${combination.MAX_TIMESTEP_SIZE}
VELOCITY_TIME_STEP_FACTOR: ${combination.VELOCITY_TIME_STEP_FACTOR}
CONSTRAINT_ITERATIONS: 4

FLUID_REST_DENSITY: 997
ETA: 1.2
EPSILON: 0.5
STIFFNESS_B: ${combination.STIFFNESS_B}
VISCOSITY_V_f: 0.05
VISCOSITY_V_b: 0.00
GRAVITY: 9.81
MAX_FLUID_PARTICLE_VELOCITY: 3.0

FluidBox:
width: 1
height: ${combination.fluid_height}
depth: 1
numberOfParticles: ${combination.numberOfParticles}
offset: 0.257 1.0 0.25

BoundaryBox:
width: 3.5
height: 2.0
depth: 1.5
samplingDistance: 0.02
`;

    const scenarioFileName = `gen_explicit_maxT${combination.MAX_TIMESTEP_SIZE}_velT${combination.VELOCITY_TIME_STEP_FACTOR}_stiff${combination.STIFFNESS_B}_fheight${combination.fluid_height}_particles${combination.numberOfParticles}.txt`;

    fs.writeFileSync(scenarioFileName, explicitScenario);
}

// Implicit Scenario
for (const combination of implicit_combinations) {
    let implicitScenario = `\
TOTAL_SECONDS_TO_SIMULATE: 4.5
MAX_TIMESTEP_SIZE: ${combination.MAX_TIMESTEP_SIZE}
VELOCITY_TIME_STEP_FACTOR: ${combination.VELOCITY_TIME_STEP_FACTOR}
CONSTRAINT_ITERATIONS: ${combination.CONSTRAINT_ITERATIONS}

FLUID_REST_DENSITY: 997
ETA: 1.2
EPSILON: 0.5
STIFFNESS_B: 1000
VISCOSITY_V_f: 0.005
VISCOSITY_V_b: 0.00
GRAVITY: 9.81
MAX_FLUID_PARTICLE_VELOCITY: 3.0

FluidBox:
width: 1
height: ${combination.fluid_height}
depth: 1
numberOfParticles: ${combination.numberOfParticles}
offset: 0.257 1.0 0.25

BoundaryBox:
width: 3.5
height: 2.0
depth: 1.5
samplingDistance: 0.02
`;

    const scenarioFileName = `gen_implicit_maxT${combination.MAX_TIMESTEP_SIZE}_velT${combination.VELOCITY_TIME_STEP_FACTOR}_iters${combination.CONSTRAINT_ITERATIONS}_fheight${combination.fluid_height}_particles${combination.numberOfParticles}.txt`;

    fs.writeFileSync(scenarioFileName, implicitScenario);
}
