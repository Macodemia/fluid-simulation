#pragma once
#include <algorithm>
#include <Eigen/Dense>
#include "scenario.h"
#include "kernel.h"
#include "data.h"

namespace learnSPH {
// -------------------------------------------------------------------------
// XSPH + Viscosity + External Accelerations
// -------------------------------------------------------------------------
inline Eigen::Vector3d getXSPHSmoothedVelocity(const Scenario& scenario, const PhysicalData& physicalData, const Neighborhood& neighborhood, const unsigned int fluidParticleIndex) {
    const Eigen::Vector3d fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    const Eigen::Vector3d fluidVel = physicalData.velocities[fluidParticleIndex];
    const double fluidDensity = physicalData.densities[fluidParticleIndex];

    Eigen::Vector3d smoothing = Eigen::Vector3d(0, 0, 0);

    neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
        const Eigen::Vector3d neighborPos = physicalData.fluidPositions[neighborIndex]; 
        const Eigen::Vector3d neighborVel = physicalData.velocities[neighborIndex];
        const double neighborDensity = physicalData.densities[neighborIndex];

        smoothing += 2.0 * physicalData.fluidParticleMass * ((neighborVel - fluidVel) / (fluidDensity + neighborDensity)) * kernel::kernelCubicFunction(fluidPos, neighborPos);
    });

    return physicalData.velocities[fluidParticleIndex] + scenario.EPSILON * smoothing;
}

inline Eigen::Vector3d getViscosityAccelerationContributionFluidToFluid(const PhysicalData& physicalData, const unsigned int fluidParticleIndex, const unsigned int neighborIndex) {
    const Eigen::Vector3d fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    const Eigen::Vector3d fluidVel = physicalData.velocities[fluidParticleIndex];
    const Eigen::Vector3d neighborPos = physicalData.fluidPositions[neighborIndex];
    const Eigen::Vector3d neighborVel = physicalData.velocities[neighborIndex];
    const double neighborDensity = physicalData.densities[neighborIndex];

    const Eigen::Vector3d fluidPosMinusNeighborPos = fluidPos - neighborPos;
    const double fluidPosMinusNeighborPosNorm = fluidPosMinusNeighborPos.norm();

    return (physicalData.fluidParticleMass / neighborDensity) * (fluidVel - neighborVel) * ( fluidPosMinusNeighborPos.transpose() * kernel::kernelGradCubicFunction(fluidPos, neighborPos, physicalData.h) / (fluidPosMinusNeighborPosNorm * fluidPosMinusNeighborPosNorm + 0.01 * physicalData.h * physicalData.h));
}

inline Eigen::Vector3d getViscosityAccelerationContributionFluidToBoundary(const PhysicalData& physicalData, const unsigned int fluidParticleIndex, const unsigned int neighborIndex) {
    const Eigen::Vector3d fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    const Eigen::Vector3d fluidVel = physicalData.velocities[fluidParticleIndex];
    const Eigen::Vector3d neighborPos = physicalData.boundaryPositions[neighborIndex];
    const double neighborVolume = physicalData.boundaryVolumes[neighborIndex];

    const Eigen::Vector3d fluidPosMinusNeighborPos = fluidPos - neighborPos;
    const double fluidPosMinusNeighborPosNorm = fluidPosMinusNeighborPos.norm();

    return neighborVolume * fluidVel * ( fluidPosMinusNeighborPos.transpose() * kernel::kernelGradCubicFunction(fluidPos, neighborPos, physicalData.h) / (fluidPosMinusNeighborPosNorm * fluidPosMinusNeighborPosNorm + 0.01 * physicalData.h * physicalData.h) );
}

inline Eigen::Vector3d getExternalAcceleration(const Scenario& scenario) {
    return Eigen::Vector3d(0, -scenario.GRAVITY, 0);
}

// -------------------------------------------------------------------------
// Weakly Compressible SPH (WCSPH)
// -------------------------------------------------------------------------
inline Eigen::Vector3d getPressureAccelerationContributionFluidToFluid(const Scenario& scenario, const PhysicalData& physicalData, const unsigned int fluidParticleIndex, const double fluidPressure, const unsigned int neighborIndex) {
    const Eigen::Vector3d fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    const double fluidDensity = physicalData.densities[fluidParticleIndex];
    const Eigen::Vector3d neighborPos = physicalData.fluidPositions[neighborIndex];
    const double neighborDensity = physicalData.densities[neighborIndex];

    const double neighborPressure = std::max(0.0, scenario.STIFFNESS_B * (neighborDensity - scenario.FLUID_REST_DENSITY));

    return physicalData.fluidParticleMass * ( (fluidPressure / (fluidDensity * fluidDensity)) + (neighborPressure / (neighborDensity * neighborDensity)) ) * kernel::kernelGradCubicFunction(fluidPos, neighborPos, physicalData.h);
}

inline Eigen::Vector3d getPressureAccelerationContributionFluidToBoundary(const Scenario& scenario, const PhysicalData& physicalData, const unsigned int fluidParticleIndex, const double fluidPressure, const unsigned int neighborIndex) {
    const Eigen::Vector3d fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    const double fluidDensity = physicalData.densities[fluidParticleIndex];
    const Eigen::Vector3d neighborPos = physicalData.boundaryPositions[neighborIndex];
    const double boundaryVolume = physicalData.boundaryVolumes[neighborIndex];

    return scenario.FLUID_REST_DENSITY * boundaryVolume * (fluidPressure / (fluidDensity * fluidDensity)) * kernel::kernelGradCubicFunction(fluidPos, neighborPos, physicalData.h);
}

// -------------------------------------------------------------------------
// Position Based Fluids (PBF)
// -------------------------------------------------------------------------
inline double computePBFSValue(const Scenario& scenario, const PhysicalData& physicalData, const Neighborhood& neighborhood, const unsigned int fluidParticleIndex) {
    const Eigen::Vector3d& fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    
    Eigen::Vector3d firstPart = Eigen::Vector3d(0, 0, 0);
    double secondPart = 0;

    // FLUID NEIGHBORS
    neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
        const Eigen::Vector3d& neighborPos = physicalData.fluidPositions[neighborIndex];
        const Eigen::Vector3d kernelGradResult = learnSPH::kernel::kernelGradCubicFunction(fluidPos, neighborPos, physicalData.h);

        firstPart += (physicalData.fluidParticleMass / scenario.FLUID_REST_DENSITY) * kernelGradResult;

        const double secondPartNorm = (-1.0 * (physicalData.fluidParticleMass / scenario.FLUID_REST_DENSITY) * kernelGradResult).norm();
        secondPart += (1.0 / physicalData.fluidParticleMass) * (secondPartNorm * secondPartNorm);
    });

    // BOUNDARY NEIGHBORS
    neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](const unsigned int neighborIndex) -> void {
        const Eigen::Vector3d& neighborPos = physicalData.boundaryPositions[neighborIndex];
        const Eigen::Vector3d kernelGradResult = learnSPH::kernel::kernelGradCubicFunction(fluidPos, neighborPos, physicalData.h);

        firstPart += physicalData.boundaryVolumes[neighborIndex] * kernelGradResult;
    });

    return (1.0 / physicalData.fluidParticleMass) * firstPart.norm() * firstPart.norm() + secondPart;
}

inline static double evaluatePBFConstraint(const Scenario& scenario, const PhysicalData& physicalData, const unsigned int fluidParticleIndex) {
    return (physicalData.densities[fluidParticleIndex] / scenario.FLUID_REST_DENSITY) - 1.0;
}

inline static double getPBFLambda(const Scenario& scenario, const PhysicalData& physicalData, const double Si, const unsigned int fluidParticleIndex) {
    const double constraint = evaluatePBFConstraint(scenario, physicalData, fluidParticleIndex);
    return constraint > 0 ? -(constraint / (Si + 1e-4)) : 0;
}

inline Eigen::Vector3d getPBFDeltaPos(const Scenario& scenario, const PhysicalData& physicalData, const Neighborhood& neighborhood, const std::vector<double>& sValues, const unsigned int fluidParticleIndex) {
    Eigen::Vector3d deltaPos = Eigen::Vector3d(0, 0, 0);
    const double lambda_i = getPBFLambda(scenario, physicalData, sValues[fluidParticleIndex], fluidParticleIndex);

    // FLUID NEIGHBORS
    neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
        const double lambda_j = getPBFLambda(scenario, physicalData, sValues[neighborIndex], neighborIndex);
        deltaPos += ((physicalData.fluidParticleMass / physicalData.fluidParticleMass) * lambda_i + lambda_j) * learnSPH::kernel::kernelGradCubicFunction(physicalData.fluidPositions[fluidParticleIndex], physicalData.fluidPositions[neighborIndex], physicalData.h);
    });
    deltaPos *= 1.0 / scenario.FLUID_REST_DENSITY;

    // BOUNDARY NEIGHBORS
    neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](const unsigned int neighborIndex) -> void {
        deltaPos += (physicalData.boundaryVolumes[neighborIndex] / physicalData.fluidParticleMass) * lambda_i * learnSPH::kernel::kernelGradCubicFunction(physicalData.fluidPositions[fluidParticleIndex], physicalData.boundaryPositions[neighborIndex], physicalData.h);
    });
    
    return deltaPos;
}

inline Eigen::Vector3d getPBFVelocityCorrection(const Eigen::Vector3d oldPosition, const Eigen::Vector3d newPosition, const double timestepSize) {
    return (newPosition - oldPosition) / timestepSize;
}

// -------------------------------------------------------------------------
// Surface Tension
// -------------------------------------------------------------------------
inline Eigen::Vector3d getCohesionForce(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const unsigned int fluidParticleIndex, const unsigned int neighborIndex) {
    if (fluidParticleIndex == neighborIndex) {
        return Eigen::Vector3d(0, 0, 0);
    }

    const Eigen::Vector3d fluidPos = physicalData.fluidPositions[fluidParticleIndex];
    const Eigen::Vector3d neighborPos = physicalData.fluidPositions[neighborIndex];
    
    const double r = (fluidPos - neighborPos).norm();
    return -scenario.COHESION * physicalData.fluidParticleMass * physicalData.fluidParticleMass * learnSPH::kernel::cohesionKernel(r, 2 * physicalData.h) * ((fluidPos - neighborPos) / (r));
}

inline void computeFluidSurfaceNormals(const learnSPH::ConfigData& configData, learnSPH::PhysicalData& physicalData, const Neighborhood& neighborhood) {
    if (physicalData.fluidSurfaceNormals.size() != physicalData.fluidPositions.size()) {
        physicalData.fluidSurfaceNormals.resize(physicalData.fluidPositions.size(), Eigen::Vector3d(0, 0, 0));
    }
    
    learnSPH::forEachParticleDo(configData.numberOfThreads, physicalData.fluidPositions, [&](const unsigned int fluidParticleIndex) -> void {
        Eigen::Vector3d normal = Eigen::Vector3d(0, 0, 0);
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
            normal += ((physicalData.fluidParticleMass) / (physicalData.densities[neighborIndex])) * learnSPH::kernel::kernelGradCubicFunction(physicalData.fluidPositions[fluidParticleIndex], physicalData.fluidPositions[neighborIndex], physicalData.h);
        });
        normal *= 2.0 * physicalData.h;

        physicalData.fluidSurfaceNormals[fluidParticleIndex] = normal;
    });
}

inline Eigen::Vector3d getCurvatureForce(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const unsigned int fluidParticleIndex, const unsigned int neighborIndex) {
    return -scenario.COHESION * physicalData.fluidParticleMass * (physicalData.fluidSurfaceNormals[fluidParticleIndex] - physicalData.fluidSurfaceNormals[neighborIndex]);
}

inline double getFinalCohesionFactor(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const unsigned int fluidParticleIndex, const unsigned int neighborIndex) {
    return (2 * scenario.FLUID_REST_DENSITY) / (physicalData.densities[fluidParticleIndex] + physicalData.densities[neighborIndex]);
}

inline Eigen::Vector3d getAdhesionForce(const Scenario& scenario, const learnSPH::PhysicalData& physicalData, const unsigned int fluidParticleIndex, const unsigned int neighborIndex) {
    const double r = (physicalData.fluidPositions[fluidParticleIndex] - physicalData.boundaryPositions[neighborIndex]).norm();

    return -scenario.ADHESION * physicalData.fluidParticleMass * physicalData.boundaryMasses[neighborIndex] * learnSPH::kernel::adhesionKernel(r, 2.0*physicalData.h) * ((physicalData.fluidPositions[fluidParticleIndex] - physicalData.boundaryPositions[neighborIndex]) / (r));
}
}
