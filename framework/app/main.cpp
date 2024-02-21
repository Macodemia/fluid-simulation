#include <iostream>
#include <thread>
#include <cmath>
#include <string>
#include "../learnSPH/common.h"
#include "../learnSPH/data.h"
#include "../learnSPH/scenario.h"
#include "../learnSPH/neighborhood.h"
#include "../learnSPH/util/vtk_writer.h"
#include "../learnSPH/physics.h"
#include "../learnSPH/particle.h"

// Reconstruction File Format:
// scenarioName
// h
// particle radius
// fluidParticles[0] -> x y z
// fluidParticles[1]
// ...
// fluidParticles[fluidParticles.size() - 1]
void saveFrameForReconstruction(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const std::vector<Eigen::Vector3d> fluidPositions, const learnSPH::ProgressData& progressData) {
    std::string reconstructionFileName = "../res/" + scenario.name + "/fluid_for_reconstruction_" + scenario.getSolverTypeString() + "_" + std::to_string(progressData.savedFrames);
    std::ofstream outFile;
    outFile.open(reconstructionFileName);

    if (!outFile.is_open()) {
        std::cout << "Error: Could not open output file for frame reconstruction. Filename: " << reconstructionFileName << std::endl;
    }

    outFile << scenario.name << std::endl;
    outFile << physicalData.h << std::endl;
    outFile << physicalData.particleRadius;
    for (const Eigen::Vector3d& particle : fluidPositions) {
        outFile << std::endl << particle.x() << " " << particle.y() << " " << particle.z();
    }

    outFile.close();
}

void saveDensityMetrics(const Scenario& scenario, const double maxDensity, const learnSPH::ProgressData& progressData) {
    std::string densityMetricsFileName = "../res/" + scenario.name + "/density_metrics_" + scenario.getSolverTypeString() + ".csv";
    std::ofstream outFile;
    outFile.open(densityMetricsFileName, std::ios_base::app);

    if (!outFile.is_open()) {
        std::cout << "Error: Could not open output file for density metrics. Filename: " << densityMetricsFileName << std::endl;
    }
    
    // Write header on first frame
    if (progressData.savedFrames == 0) {
        outFile << "Difference to Rest Density" << std::endl;
    }

    outFile << maxDensity - scenario.FLUID_REST_DENSITY << std::endl;

    outFile.close();
}

void saveFluidVTK(const Scenario& scenario, const std::vector<Eigen::Vector3d>& fluidPositions, const std::vector<double>& fluidDensities, const learnSPH::ProgressData& progressData) {
    const std::string filenameFluid = "../res/" + scenario.name + "/fluid_" + scenario.getSolverTypeString() + "_" + std::to_string(progressData.savedFrames) + ".vtk";
    const std::string dirName = std::filesystem::path(filenameFluid).parent_path();
    std::filesystem::create_directories(dirName);
    std::vector<Eigen::Vector3d> vectorDataFluidParticles(fluidPositions.size(), Eigen::Vector3d(0,0,0));
    learnSPH::save_particles_to_vtk(filenameFluid, fluidPositions, fluidDensities, vectorDataFluidParticles);
}

void saveBoundaryVTK(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const learnSPH::ProgressData& progressData) {
    for (unsigned int groupIndex = 0; groupIndex < physicalData.boundaryGroups.size(); groupIndex++) {
        const learnSPH::BoundaryParticleGroup& boundaryGroup = physicalData.boundaryGroups[groupIndex];
        std::vector<Eigen::Vector3d> groupParticles;
        groupParticles.reserve(boundaryGroup.lastIndex - boundaryGroup.firstIndex + 1);
        for (unsigned int i = boundaryGroup.firstIndex; i <= boundaryGroup.lastIndex; i++) {
            groupParticles.push_back(physicalData.boundaryPositions[i]); 
        }

        const std::string filenameBoundary = "../res/" + scenario.name + "/boundary_" + scenario.getSolverTypeString() + "_group" + std::to_string(groupIndex) + "_" + std::to_string(progressData.savedFrames) + ".vtk";
        const std::string dirName = std::filesystem::path(filenameBoundary).parent_path();
        std::filesystem::create_directories(dirName);
        std::vector<Eigen::Vector3d> vectorDataBoundaryParticles(groupParticles.size(), Eigen::Vector3d(0,0,0));
        std::vector<double> boundaryParticleDensity(groupParticles.size(), 0);
        learnSPH::save_particles_to_vtk(filenameBoundary, groupParticles, boundaryParticleDensity, vectorDataBoundaryParticles);
    }
}

void saveVTK(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const std::vector<Eigen::Vector3d>& fluidPositions, const std::vector<double>& fluidDensities, const learnSPH::ProgressData& progressData) {
    saveFluidVTK(scenario, fluidPositions, fluidDensities, progressData);
    if (!physicalData.boundaryPositions.empty() && progressData.savedFrames == 0) {
        saveBoundaryVTK(scenario, physicalData, progressData);
    }
}

inline double getNextTimestepSize(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const learnSPH::ProgressData& progressData) {
    double maxVelocityValue = -std::numeric_limits<double>::max();
    for (auto& velocity : physicalData.velocities) {
        maxVelocityValue = std::max(maxVelocityValue, velocity.norm());
    }

    // upper-bounded by maxTimeStepSize
    double timestepSize = std::min(scenario.VELOCITY_TIME_STEP_FACTOR * (physicalData.particleRadius / maxVelocityValue), scenario.MAX_TIMESTEP_SIZE);

    return timestepSize;
}

// First argument is either "explicit" or "implicit" -> Specifies solver
// Second (optional) argument is a scenario file
int main(int argc, char** argv)
{
    Scenario scenario;
    if (argc < 2) {
        std::cout << "ERROR: Please specify the solver!" << std::endl;
        std::cout << "Usage (script): ./run(_release).sh implicit|explicit scenarios/example_name.txt" << std::endl;
        std::cout << "Usage (direct): ./build/xxx/app/learnSPH implicit|explicit scenarios/example_name.txt [threadCount]" << std::endl;
        return 1;
    } else if (argc < 3) {
        std::cout << "ERROR: No scenario file specified!" << std::endl;
        std::cout << "Usage (script): ./run(_release).sh implicit|explicit scenarios/example_name.txt" << std::endl;
        std::cout << "Usage (direct): ./build/xxx/app/learnSPH implicit|explicit scenarios/example_name.txt [threadCount]" << std::endl;
        return 1;
    }
    
    // The solver must be set before reading a scenario file!
    std::string solver(argv[1]);
    if (solver == "explicit") {
        scenario.solverType = SOLVER_TYPE::EXPLICIT;
    } else if (solver == "implicit") {
        scenario.solverType = SOLVER_TYPE::IMPLICIT;
    } else {
        std::cout << "ERROR: Solver must either be explicit or implicit! Provided solver: " << solver << std::endl;
        std::cout << "Usage (script): ./run(_release).sh implicit|explicit scenarios/example_name.txt [threadCount]" << std::endl;
        std::cout << "Usage (direct): ./build/xxx/app/learnSPH implicit|explicit scenarios/example_name.txt [threadCount]" << std::endl;
        return 1;
    }

    std::string scenarioFile(argv[2]);
    if (std::filesystem::exists(scenarioFile) && std::filesystem::is_regular_file(scenarioFile)) {
        scenario.readFromFile(scenarioFile);
        if (scenario.fluidBoxes.empty()) {
            std::cout << "ERROR: Specified scenario (" << scenarioFile << ") does not specify a fluid box!" << std::endl;
            return 1;
        }
    } else {
        std::cout << "ERROR: Specified scenario file (" << scenarioFile << ") does not exist!" << std::endl;
        return 1;
    }

    std::cout << std::endl;
    scenario.print();

    learnSPH::ConfigData configData;
    // Read number of threads if given
    if (argc >= 4) {
        configData.numberOfThreads = std::stoi(argv[3]);   
    } else {
        configData.numberOfThreads = std::thread::hardware_concurrency();
    }
    configData.print();

    learnSPH::PhysicalData physicalData;
    learnSPH::ProgressData progressData;

    // -------------------------------------------------------------------------
    // Sample Fluid Particles
    //     TODO: Currently only one fluid box is allowed
    // -------------------------------------------------------------------------
    const FluidBoxParameter& fluidBox = scenario.fluidBoxes[0];
    learnSPH::particle::sampleParticles(physicalData, fluidBox.width, fluidBox.height, fluidBox.depth, scenario.FLUID_REST_DENSITY, fluidBox.numberOfParticles, fluidBox.offset, fluidBox.rotation);

    // -------------------------------------------------------------------------
    // Sample Boundary Particles
    // -------------------------------------------------------------------------
    // From boundary boxes
    for (const BoundaryBoxParameter& boundaryBox : scenario.boundaryBoxes) {
        const std::vector<Eigen::Vector3d> sampledBoundaryParticles = learnSPH::particle::sampleBoundaryBox(boundaryBox.width, boundaryBox.height, boundaryBox.depth, boundaryBox.samplingDistance, boundaryBox.offset);
        if (sampledBoundaryParticles.empty()) {
            continue;
        }
        physicalData.boundaryGroups.push_back(learnSPH::BoundaryParticleGroup(physicalData.boundaryPositions.size(), physicalData.boundaryPositions.size() + sampledBoundaryParticles.size() - 1));
        
        physicalData.boundaryPositions.insert(physicalData.boundaryPositions.end(), sampledBoundaryParticles.begin(), sampledBoundaryParticles.end());
    }

    // From boundary meshes
    for (const BoundaryMeshParameter& boundaryMesh : scenario.boundaryMeshes) {
        const std::vector<Eigen::Vector3d> sampledBoundaryParticles = learnSPH::particle::sampleObjFile(boundaryMesh.file, boundaryMesh.samplingDistance, boundaryMesh.offset, boundaryMesh.scale);
        if (sampledBoundaryParticles.empty()) {
            continue;
        }
        physicalData.boundaryGroups.push_back(learnSPH::BoundaryParticleGroup(physicalData.boundaryPositions.size(), physicalData.boundaryPositions.size() + sampledBoundaryParticles.size() - 1));

        physicalData.boundaryPositions.insert(physicalData.boundaryPositions.end(), sampledBoundaryParticles.begin(), sampledBoundaryParticles.end());
    }

    physicalData.boundaryMasses.resize(physicalData.boundaryPositions.size(), 0);
    physicalData.boundaryVolumes.resize(physicalData.boundaryPositions.size(), 0);

    // -------------------------------------------------------------------------
    // Various Initializations
    //     Neighborhood creation
    //     h calculation
    //     FLUID_REST_DENSITY correction
    // -------------------------------------------------------------------------
    physicalData.computeH(scenario.FLUID_REST_DENSITY, scenario.ETA);

    // NOTE: This is just a temporary neighborhood, the real neighborhood with the corrected h will be created after rest density correction
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);

    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, scenario.FLUID_REST_DENSITY, configData.numberOfThreads);
    
    // FLUID_REST_DENSITY Correction to prohibit initial "explosion" -> Search max estimated density -> set FLUID_REST_DENSITY to it
    {
        double maxEstimatedDensity = 0;
        for (const double density : physicalData.densities) {
            maxEstimatedDensity = std::max(maxEstimatedDensity, density);
        }
        scenario.FLUID_REST_DENSITY = maxEstimatedDensity * 1.05; // TODO: What is a good factor here?
    }

    physicalData.computeH(scenario.FLUID_REST_DENSITY, scenario.ETA);
    neighborhood = learnSPH::Neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    neighborhood.update();

    // -------------------------------------------------------------------------
    // Simulation Loop
    // -------------------------------------------------------------------------
    std::vector<Eigen::Vector3d> accelerations(physicalData.fluidPositions.size(), Eigen::Vector3d(0, 0, 0));
    std::vector<Eigen::Vector3d> previousFluidPositions = physicalData.fluidPositions;
    std::vector<double> previousDensities = physicalData.densities;
    std::vector<double> sValues(physicalData.fluidPositions.size(), 0);
   
    // Save start frame
    saveVTK(scenario, physicalData, physicalData.fluidPositions, physicalData.densities, progressData);
    saveFrameForReconstruction(scenario, physicalData, physicalData.fluidPositions, progressData);
    double maxDensity = 0;
    // NOTE: Using only one thread, as otherwise there would be data races
    learnSPH::forEachParticleDo(1, physicalData.fluidPositions, [&](const unsigned int i) -> void {
        maxDensity = std::max(maxDensity, physicalData.densities[i]);
    });
    saveDensityMetrics(scenario, maxDensity, progressData);
    progressData.savedFrames++;
    
    physicalData.print(scenario);
    std::cout << std::endl;
    const auto beginOfSimulation = std::chrono::high_resolution_clock::now();
    double averageFrameTime = 0;
    while (progressData.simulatedTime < scenario.TOTAL_SECONDS_TO_SIMULATE) {
        // -------------------------------------------------------------------------
        // Loop Preperation
        // -------------------------------------------------------------------------
        const auto beginOfFrame = std::chrono::high_resolution_clock::now();
        double timestepSize = getNextTimestepSize(scenario, physicalData, progressData);
        
        if (scenario.solverType == SOLVER_TYPE::EXPLICIT) {
            neighborhood.update();
        }
        physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, scenario.FLUID_REST_DENSITY, configData.numberOfThreads);
        
        // -------------------------------------------------------------------------
        // Loop Physics
        // -------------------------------------------------------------------------
        
        // Precompute Fluid Surface Normals
        learnSPH::computeFluidSurfaceNormals(configData, physicalData, neighborhood);

        // Particle Loop - Compute Accelerations
        learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
            // External Acceleration
            accelerations[fluidParticleIndex] = learnSPH::getExternalAcceleration(scenario); // IMPORTANT: Use = (equals) not += (plus equals)! We do not want to accumulate accelerations over multiple frames
          
            double fluidPressure;
            if (scenario.solverType == SOLVER_TYPE::EXPLICIT) {
                fluidPressure = std::max(0.0, scenario.STIFFNESS_B * (physicalData.densities[fluidParticleIndex] - scenario.FLUID_REST_DENSITY));
            }

            // FLUID NEIGHBORS
            neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
                // Viscosity Acceleration
                accelerations[fluidParticleIndex] += 2.0 * scenario.VISCOSITY_V_f * learnSPH::getViscosityAccelerationContributionFluidToFluid(physicalData, fluidParticleIndex, neighborIndex);

                // Cohesion + Curvature Factor
                const double cohesionAndCurvatureFactor = learnSPH::getFinalCohesionFactor(scenario, physicalData, fluidParticleIndex, neighborIndex);
                // Cohesion Force
                accelerations[fluidParticleIndex] += cohesionAndCurvatureFactor * learnSPH::getCohesionForce(scenario, physicalData, fluidParticleIndex, neighborIndex) / physicalData.fluidParticleMass;
                // Curvature Force
                accelerations[fluidParticleIndex] += cohesionAndCurvatureFactor * learnSPH::getCurvatureForce(scenario, physicalData, fluidParticleIndex, neighborIndex) / physicalData.fluidParticleMass;

                if (scenario.solverType == SOLVER_TYPE::EXPLICIT) {
                    accelerations[fluidParticleIndex] -= learnSPH::getPressureAccelerationContributionFluidToFluid(scenario, physicalData, fluidParticleIndex, fluidPressure, neighborIndex);
                }
            });

            // BOUNDARY NEIGHBORS
            neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](const unsigned int neighborIndex) -> void {
                // Viscosity Acceleration
                accelerations[fluidParticleIndex] += 2.0 * scenario.VISCOSITY_V_b * learnSPH::getViscosityAccelerationContributionFluidToBoundary(physicalData, fluidParticleIndex, neighborIndex);

                // Adhesion Force
                accelerations[fluidParticleIndex] += learnSPH::getAdhesionForce(scenario, physicalData, fluidParticleIndex, neighborIndex) / physicalData.fluidParticleMass;

                if (scenario.solverType == SOLVER_TYPE::EXPLICIT) {
                    accelerations[fluidParticleIndex] -= learnSPH::getPressureAccelerationContributionFluidToBoundary(scenario, physicalData, fluidParticleIndex, fluidPressure, neighborIndex);
                }
            });
        });

        // Particle Loop - Compute Velocities
        learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
            physicalData.velocities[fluidParticleIndex] += timestepSize * accelerations[fluidParticleIndex];
            // Cap to max velocity
            const double velocityMagnitude = physicalData.velocities[fluidParticleIndex].norm();
            if (velocityMagnitude > scenario.MAX_FLUID_PARTICLE_VELOCITY) {
                physicalData.velocities[fluidParticleIndex].normalize();
                physicalData.velocities[fluidParticleIndex] *= scenario.MAX_FLUID_PARTICLE_VELOCITY;
            }
        });

        // Particle Loop - Smoothed Velocities
        std::vector<Eigen::Vector3d> smoothedVelocities(physicalData.velocities.size());
        learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
            smoothedVelocities[fluidParticleIndex] = learnSPH::getXSPHSmoothedVelocity(scenario, physicalData, neighborhood, fluidParticleIndex);
        });

        // Particle Loop - Compute Positions
        learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
            physicalData.fluidPositions[fluidParticleIndex] += timestepSize * smoothedVelocities[fluidParticleIndex];
        });

        if (scenario.solverType == SOLVER_TYPE::IMPLICIT) {
            neighborhood.update();

            // PBF Iterations
            for (unsigned int iter = 0; iter < scenario.CONSTRAINT_ITERATIONS; iter++) {
                physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, scenario.FLUID_REST_DENSITY, configData.numberOfThreads);

                // Compute S_i Values
                learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
                    sValues[fluidParticleIndex] = learnSPH::computePBFSValue(scenario, physicalData, neighborhood, fluidParticleIndex);
                });

                // Compute delta positions (delta x)
                std::vector<Eigen::Vector3d> deltaPositions(physicalData.fluidPositions.size());
                learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
                    deltaPositions[fluidParticleIndex] = learnSPH::getPBFDeltaPos(scenario, physicalData, neighborhood, sValues, fluidParticleIndex);
                });

                // Update Particle Position
                learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
                    physicalData.fluidPositions[fluidParticleIndex] += deltaPositions[fluidParticleIndex];
                });
            }

            // Correct Velocities
            learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
                // NOTE: The positions are already correct as we used physicalData.fluidPositions above. We can do that because we have the previousFluidPositions anyway (Because we need them for interpolation)
                physicalData.velocities[fluidParticleIndex] = learnSPH::getPBFVelocityCorrection(previousFluidPositions[fluidParticleIndex], physicalData.fluidPositions[fluidParticleIndex], timestepSize);
            });
        }

        // -------------------------------------------------------------------------
        // Loop Completion
        // -------------------------------------------------------------------------

        // Save frame if necessary
        if (progressData.timeSinceLastSave >= scenario.SAVE_FRAME_INTERVAL) {
            // Interpolation using a+t(b−a) - Calculate position at the correct time by interpolating between last positions and new positions
            // NOTE: Reuse previousFluidPositions to avoid allocation as it will be overwritten after frame saving anyways
            // NOTE: Reuse previousDensities to avoid allocation as it will be overwritten after frame saving anyways
            const double t = 1 - ((progressData.timeSinceLastSave - scenario.SAVE_FRAME_INTERVAL) / timestepSize);
            std::atomic<double> maxDensity = 0;
            learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int i) -> void {
                previousFluidPositions[i] = previousFluidPositions[i] + t * (physicalData.fluidPositions[i] - previousFluidPositions[i]);
                previousDensities[i] = previousDensities[i] + t * (physicalData.densities[i] - previousDensities[i]);

                maxDensity = std::max(previousDensities[i], maxDensity.load());
            });

            saveVTK(scenario, physicalData, previousFluidPositions, previousDensities, progressData);
            saveDensityMetrics(scenario, maxDensity, progressData);
            saveFrameForReconstruction(scenario, physicalData, previousFluidPositions, progressData);

            progressData.savedFrames++;
            progressData.timeSinceLastSave = 0;
        }

        previousFluidPositions = physicalData.fluidPositions; // NOTE: This must be done after the frame got saved (Save Frames uses the previousFluidPositions)
        previousDensities = physicalData.densities; // NOTE: This must be done after the frame got saved (Save Frames uses the previousDensities)

        // Update timing
        progressData.frames++;
        progressData.simulatedTime += timestepSize;
        progressData.timeSinceLastSave += timestepSize;

        const auto endOfFrame = std::chrono::high_resolution_clock::now();
        const double frameDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endOfFrame-beginOfFrame).count();
        const int runtime = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - beginOfSimulation).count();
        averageFrameTime += frameDuration;
        const double currentAverageFrameTime = averageFrameTime / progressData.frames;
        const double estimatedRuntimeRemaining = (scenario.TOTAL_SECONDS_TO_SIMULATE - progressData.simulatedTime);
        
        // \e[2K clears the current line and \r sets the cursor to the beginning of the last printed line
        // The std::flush is needed to see the line (normally std::endl flushes the stdout buffer)
        std::cout << "\e[2K" << "\r";
        std::string frameDurationString = std::to_string(frameDuration).substr(0, std::to_string(frameDuration).find('.'));
        std::string outputString = "Remaining Time: " + std::to_string(scenario.TOTAL_SECONDS_TO_SIMULATE - progressData.simulatedTime)
            + " | Saved Frames: " + std::to_string(progressData.savedFrames)
            + " | Timestep Size: " + std::to_string(timestepSize) + "ms"
            + " | Frame Computation Time: " + frameDurationString + "ms"
            + " | Runtime: " + std::to_string(runtime) + "s";
        printf("| %-*s |", 50, outputString.c_str());
        std::cout << std::flush;
    }

    // Print Metrics
    averageFrameTime /= progressData.frames;
    const auto endOfSimulation = std::chrono::high_resolution_clock::now();
    const double completeDurationInMs = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSimulation-beginOfSimulation).count();

    std::cout << std::endl << std::endl;
    std::cout << "--------------------------" << std::endl;
    std::cout << "Finished:" << std::endl;
    std::cout << "\tCompleted in " << completeDurationInMs / 1000.0 << " seconds!" << std::endl;
    std::cout << "\tAverage Frame Computation Time: " << averageFrameTime << "ms" << std::endl;
    std::cout << "--------------------------" << std::endl;
    return 0;
}
