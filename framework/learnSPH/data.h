#pragma once
#include <iostream>
#include "common.h"
#include "../learnSPH/kernel.h"
#include "neighborhood.h"
#include "scenario.h"
#include <Eigen/Dense>

namespace learnSPH {

struct ConfigData {
    unsigned int numberOfThreads = 1;

    void print() const {
        std::cout << "--------------------------" << std::endl;
        std::cout << "Config:" << std::endl;
        std::cout << "\tNumber of Threads: " << numberOfThreads << std::endl;
        std::cout << "--------------------------" << std::endl;
    }
};

struct BoundaryParticleGroup {
    unsigned int firstIndex;
    unsigned int lastIndex;

    BoundaryParticleGroup(unsigned int startIndex, unsigned int endIndex) {
        this->firstIndex = startIndex;
        this->lastIndex = endIndex;
    }
};

struct PhysicalData {
    double particleRadius;
    double h;
    double fluidParticleMass;
    std::vector<Eigen::Vector3d> fluidPositions;
    std::vector<Eigen::Vector3d> velocities;
    std::vector<Eigen::Vector3d> fluidSurfaceNormals;
    std::vector<double> densities;
    
    // Boundary Particles
    std::vector<BoundaryParticleGroup> boundaryGroups;
    std::vector<Eigen::Vector3d> boundaryPositions;
    std::vector<double> boundaryVolumes;
    std::vector<double> boundaryMasses;
    bool boundaryVolumesNeedUpdate = true;
    
    void updateParticleDensitiesAndBoundaryVolumes(const learnSPH::Neighborhood& neighborhood, const double fluidRestDensity, const unsigned int numberOfThreads) {
        const double boundaryRestDensity = fluidRestDensity; // TODO: What is correct here?

        // Boundary particle volumes only need to be updated when they moved
        if (boundaryVolumesNeedUpdate && !boundaryVolumes.empty()) {
            learnSPH::forEachParticleDo(numberOfThreads, boundaryPositions, [&](unsigned int boundaryParticleIndex) -> void {
                // compute the volume
                double kernelValues = 0;
                neighborhood.forEachNeighborDo(neighborhood.boundaryPointSetId, boundaryParticleIndex, neighborhood.boundaryPointSetId, [&](const unsigned int neighborIndex) -> void {
                    kernelValues += learnSPH::kernel::kernelCubicFunction(boundaryPositions[boundaryParticleIndex], boundaryPositions[neighborIndex]);
                });

                boundaryVolumes[boundaryParticleIndex] = 1.0 / kernelValues;
                boundaryMasses[boundaryParticleIndex] = boundaryRestDensity * boundaryVolumes[boundaryParticleIndex];
            });

            boundaryVolumesNeedUpdate = false;
        }

        // Compute fluid density
        const CompactNSearch::PointSet& fluidPointSet = neighborhood.nsearch.point_set(neighborhood.fluidPointSetId);
        const CompactNSearch::PointSet& boundaryPointSet = neighborhood.nsearch.point_set(neighborhood.boundaryPointSetId);

        // Compute density for fluid particles
        learnSPH::forEachParticleDo(numberOfThreads, fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
            densities[fluidParticleIndex] = 0; // IMPORTANT: We need to reset the density value here as it is not accumulated over time/frames

            // Neighbor is a fluid particle
            neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
                densities[fluidParticleIndex] += fluidParticleMass * learnSPH::kernel::kernelCubicFunction(fluidPositions[fluidParticleIndex], fluidPositions[neighborIndex]);
            });

            // Neighbor is a boundary particle
            if (!boundaryPositions.empty()) {
                neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](const unsigned int neighborIndex) -> void {
                    densities[fluidParticleIndex] += boundaryRestDensity * boundaryVolumes[neighborIndex] * learnSPH::kernel::kernelCubicFunction(fluidPositions[fluidParticleIndex], boundaryPositions[neighborIndex]);
                });
            }
        });
    }
    
    void computeH(const double fluidRestDensity, const double eta) {
        const double representativeDiameter = std::pow(fluidParticleMass / fluidRestDensity, 1.0/3.0);
        h = eta * representativeDiameter;
        learnSPH::kernel::computeLookupTables(h);
    }

    void print(const Scenario& scenario) {
        std::cout << "--------------------------" << std::endl;
        std::cout << "Physical Data:" << std::endl;
        std::cout << "\tFluid Particle Radius: " << particleRadius << std::endl;
        std::cout << "\tFluid Particle Mass: " << fluidParticleMass << std::endl;
        std::cout << "\tFluid Particle Count: " << fluidPositions.size() << std::endl;
        std::cout << "\tCorrected Fluid Rest Density: " << scenario.FLUID_REST_DENSITY << std::endl;
        std::cout << "\tSmoothing Length h: " << h << std::endl;
        std::cout << "\tBoundary Particle Count: " << boundaryPositions.size() << std::endl;
        std::cout << "--------------------------" << std::endl;
    }
};

struct ProgressData {
    uint64 frames = 0;
    double simulatedTime = 0;
    uint64 savedFrames = 0;
    double timeSinceLastSave = 0;
};

} // namespace learnSPH
