#include "catch.hpp"
#include "../learnSPH/scenario.h"
#include "../learnSPH/data.h"
#include "../learnSPH/neighborhood.h"
#include "../learnSPH/physics.h"
#include <random>
#include <algorithm>

static std::vector<Eigen::Vector3d> generateRandomVectors3d(int count) {
    const int closeVectorsCount = (int) ((float) count * 0.1f);
    assert(count - closeVectorsCount > 0);

    std::uniform_real_distribution<double> distribution(-5.0, 5.0);
    std::random_device r;
    std::default_random_engine randomEngine(r());

    std::vector<Eigen::Vector3d> result;

    for (int i = 0; i < count - closeVectorsCount; i++) {

        const double x = distribution(randomEngine);
        const double y = distribution(randomEngine);
        const double z = distribution(randomEngine);

        result.emplace_back(x, y, z);
    }

    distribution = std::uniform_real_distribution<double>(-0.5, 0.5);
    for (int i = 0; i < closeVectorsCount; i++) {
        const double x = distribution(randomEngine);
        const double y = distribution(randomEngine);
        const double z = distribution(randomEngine);

        result.emplace_back(x, y, z);
    }

    return result;
}

static std::vector<double> generateRandomDoubles(int count, double min, double max) {
    std::uniform_real_distribution<double> distribution(min, max);
    std::random_device r;
    std::default_random_engine randomEngine(r());

    std::vector<double> result;
    for (int i = 0; i < count; i++) {
        result.push_back(distribution(randomEngine));
    }

    return result;
}

static double generateRandomDouble(double min, double max) {
    std::uniform_real_distribution<double> distribution(min, max);
    std::random_device r;
    std::default_random_engine randomEngine(r());
    return distribution(randomEngine);
}

static learnSPH::PhysicalData generateRandomPhysicalData(int particleCount) {
    learnSPH::PhysicalData physicalData;
    physicalData.fluidPositions = generateRandomVectors3d(particleCount);
    physicalData.velocities = generateRandomVectors3d(particleCount);
    physicalData.densities = generateRandomDoubles(particleCount, 1, 2000);
    physicalData.fluidSurfaceNormals = generateRandomVectors3d(particleCount);
    physicalData.fluidParticleMass = generateRandomDouble(0.0001, 2.0);
    physicalData.particleRadius = generateRandomDouble(0.0001, 1.0);
    physicalData.boundaryPositions = generateRandomVectors3d(particleCount);
    physicalData.boundaryVolumes = generateRandomDoubles(particleCount, 0.0001, 10.0);
    physicalData.boundaryMasses = generateRandomDoubles(particleCount, 0.001, 10.0);
    physicalData.computeH(1000, 1.2);
    return physicalData;
}

static double distanceNorm(Eigen::Vector3d x1, Eigen::Vector3d x2) {
    return (x1-x2).norm();
}

static double s(double q) {
    return ((3.0)/(2.0*(M_PI)))*
    ((0.0 <= q && q < 1.0) ? ((2.0)/(3.0))-std::pow(q, 2)+((1.0)/(2.0))*std::pow(q, 3) : ((1.0 <= q && q < 2.0) ? ((1.0)/(6.0))*(std::pow(2-q, 3)) : 0));
}
static double W(Eigen::Vector3d x1, Eigen::Vector3d x2, double h) {
    return ((1.0)/(std::pow(h, 3)))*s((((distanceNorm(x1,x2)))/(h)));
}
static double sgrad(double q) {
    return ((0.0 <= q && q < 1.0) ? (((3.0)/(2.0*(M_PI)))) * (-2.0*q+((3.0)/(2.0))*std::pow(q, 2)) : ((1.0 <= q && q < 2.0) ? (((3.0)/(2.0*(M_PI)))) * (((-1.0)/(2.0))*(2.0-q)*(2.0-q)) : 0.0));
}
static Eigen::Vector3d Wgrad(Eigen::Vector3d x1, Eigen::Vector3d x2, double h) {
    return ((x1 == x2) ? Eigen::Vector3d(0,0,0) : Eigen::Vector3d(
            ((((1)/(std::pow(h, 4)))*sgrad((((distanceNorm(x1,x2)))/(h))))/((distanceNorm(x1,x2)))) * (x1.x() - x2.x()),
            ((((1)/(std::pow(h, 4)))*sgrad((((distanceNorm(x1,x2)))/(h))))/((distanceNorm(x1,x2)))) * (x1.y() - x2.y()),
            ((((1)/(std::pow(h, 4)))*sgrad((((distanceNorm(x1,x2)))/(h))))/((distanceNorm(x1,x2)))) * (x1.z() - x2.z())
        ));
}
static Eigen::Vector3d xsph(double m, Eigen::Vector3d x1, Eigen::Vector3d x2, Eigen::Vector3d v1, Eigen::Vector3d v2, double p1, double p2, double h) {
    return 2.0*m*((v2-v1)/(p1+p2))*W(x1,x2,h);
}
static double P_i(double B, double pi, double p0) {
    return std::max(0.0, B*(pi-p0));
}
static Eigen::Vector3d apfluidtofluid(double m, Eigen::Vector3d x1, Eigen::Vector3d x2, double pressure1, double pressure2, double density1, double density2, double h) {
    return m*(((pressure1)/(std::pow(density1, 2)))+((pressure2)/(std::pow(density2, 2)))) * Wgrad(x1, x2, h);
}
static Eigen::Vector3d apfluidtoboundary(double restDens, Eigen::Vector3d x1, Eigen::Vector3d x2, double volume, double pressure1, double density1, double h) {
    return restDens*volume*(((pressure1)/(std::pow(density1, 2)))) * Wgrad(x1, x2, h);
}
static Eigen::Vector3d aviscosfluidtofluid(Eigen::Vector3d x1, Eigen::Vector3d x2, Eigen::Vector3d v1, Eigen::Vector3d v2, double m, double density2, double h) {
    return ((m)/(density2))*(v1-v2)*((((x1-x2).transpose()*Wgrad(x1,x2,h))/((distanceNorm(x1,x2))*(distanceNorm(x1,x2))+0.01*std::pow(h, 2))));
}
static Eigen::Vector3d aviscosfluidtoboundary(Eigen::Vector3d x1, Eigen::Vector3d x2, Eigen::Vector3d v1, double volume, double h) {
    return volume*v1*((((x1-x2).transpose()*Wgrad(x1,x2,h))/((distanceNorm(x1,x2))*(distanceNorm(x1,x2))+0.01*std::pow(h, 2))));
}
static double Wcoh(double r, double c) {
    return ((32)/((M_PI)*std::pow(c, 9)))*
    ((0 <= r && r <= ((c)/(2))) ? 2*std::pow(c-r, 3)*std::pow(r, 3)-((std::pow(c, 6))/(64)) : ((((c)/(2)) < r && r <= c) ? std::pow(c-r, 3)*std::pow(r, 3) : 0));
}
static double Wadh(double r, double c) {
    return ((0.007)/(std::pow(c, 3.25))*(((c)/(2) <= r && r <= c) ? std::pow(((-4*std::pow(r, 2))/(c))+6*r-2*c, 1.0/4.0) : 0));
}
static Eigen::Vector3d fcohesion(double cohFac, double m, Eigen::Vector3d x1, Eigen::Vector3d x2, double c) {
    if (x1 == x2) {
        return Eigen::Vector3d(0,0,0);
    }
    return -cohFac*m*m*Wcoh((distanceNorm(x1,x2)),c)*((x1-x2)/((distanceNorm(x1,x2))));
}
static Eigen::Vector3d ni(double c, double m, double density, Eigen::Vector3d x1, Eigen::Vector3d x2, double h) {
    return c*((m)/(density))*Wgrad(x1,x2,h);
}
static Eigen::Vector3d fcurvature(double cohFac, double m, Eigen::Vector3d n1, Eigen::Vector3d n2) {
    return -cohFac*m*(n1-n2);
}
static Eigen::Vector3d fadhesion(double adhFac, double mFluid, double mBound, Eigen::Vector3d x1, Eigen::Vector3d x2, double c) {
    return -adhFac*mFluid*mBound*Wadh((distanceNorm(x1,x2)),c)*((x1-x2)/((distanceNorm(x1,x2))));
}

TEST_CASE("XSPH Smoothed Velocity", "[Physics]") {
    //Scenario scenario;
    //learnSPH::PhysicalData physicalData = generateRandomPhysicalData(1000);
    // TODO
}

TEST_CASE("Viscosity Acceleration Fluid To Fluid Contribution", "[Physics]") {
    Scenario scenario;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
            const Eigen::Vector3d actualResult = learnSPH::getViscosityAccelerationContributionFluidToFluid(physicalData, fluidParticleIndex, neighborIndex);
            const Eigen::Vector3d expectedResult = aviscosfluidtofluid(physicalData.fluidPositions[fluidParticleIndex], physicalData.fluidPositions[neighborIndex],
                                                                       physicalData.velocities[fluidParticleIndex], physicalData.velocities[neighborIndex],
                                                                       physicalData.fluidParticleMass, physicalData.densities[neighborIndex], physicalData.h);

            //std::cout << "ACTUAL: " << actualResult << std::endl;
            //std::cout << "EXPECTED: " << expectedResult << std::endl;
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

TEST_CASE("Viscosity Acceleration Fluid To Boundary Contribution", "[Physics]") {
    Scenario scenario;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](unsigned int neighborIndex) -> void {
            

            const Eigen::Vector3d actualResult = learnSPH::getViscosityAccelerationContributionFluidToBoundary(physicalData, fluidParticleIndex, neighborIndex);
            const Eigen::Vector3d expectedResult = aviscosfluidtoboundary(physicalData.fluidPositions[fluidParticleIndex], physicalData.boundaryPositions[neighborIndex], physicalData.velocities[fluidParticleIndex], physicalData.boundaryVolumes[neighborIndex], physicalData.h); 

            //std::cout << "ACTUAL: " << actualResult << std::endl;
            //std::cout << "EXPECTED: " << expectedResult << std::endl;
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

/*
Eigen::Vector3d getXSPHSmoothedVelocity(const Scenario& scenario, const PhysicalData& physicalData, const Neighborhood& neighborhood, const unsigned int fluidParticleIndex)
*/
TEST_CASE("WCSPH Pressure Acceleration Fluid to Fluid Contribution", "[Physics]") {
    Scenario scenario;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(100000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
            
            const double actualFluidPressure = std::max(0.0, scenario.STIFFNESS_B * (physicalData.densities[fluidParticleIndex] - scenario.FLUID_REST_DENSITY));
            const Eigen::Vector3d actualResult = learnSPH::getPressureAccelerationContributionFluidToFluid(scenario, physicalData, fluidParticleIndex, actualFluidPressure, neighborIndex);
            const double pressure1 = P_i(scenario.STIFFNESS_B, physicalData.densities[fluidParticleIndex], scenario.FLUID_REST_DENSITY);
            const double pressure2 = P_i(scenario.STIFFNESS_B, physicalData.densities[neighborIndex], scenario.FLUID_REST_DENSITY);
            const Eigen::Vector3d expectedResult = apfluidtofluid(physicalData.fluidParticleMass, physicalData.fluidPositions[fluidParticleIndex], physicalData.fluidPositions[neighborIndex], pressure1, pressure2, physicalData.densities[fluidParticleIndex], physicalData.densities[neighborIndex], physicalData.h);

            if (actualResult.norm() > 0 || expectedResult.norm() > 0) {
                //std::cout << "ACTUAL: " << actualResult << std::endl;
                //std::cout << "EXPECTED: " << expectedResult << std::endl;
            }
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

TEST_CASE("WCSPH Pressure Acceleration Fluid to Boundary Contribution", "[Physics]") {
    Scenario scenario;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(100000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](unsigned int neighborIndex) -> void {
            
            const double actualFluidPressure = std::max(0.0, scenario.STIFFNESS_B * (physicalData.densities[fluidParticleIndex] - scenario.FLUID_REST_DENSITY));
            const Eigen::Vector3d actualResult = learnSPH::getPressureAccelerationContributionFluidToBoundary(scenario, physicalData, fluidParticleIndex, actualFluidPressure, neighborIndex); 
            const double pressure1 = P_i(scenario.STIFFNESS_B, physicalData.densities[fluidParticleIndex], scenario.FLUID_REST_DENSITY);
            const Eigen::Vector3d expectedResult = apfluidtoboundary(scenario.FLUID_REST_DENSITY, physicalData.fluidPositions[fluidParticleIndex], physicalData.boundaryPositions[neighborIndex], physicalData.boundaryVolumes[neighborIndex], pressure1, physicalData.densities[fluidParticleIndex], physicalData.h);

            if (actualResult.norm() > 0 || expectedResult.norm() > 0) {
                //std::cout << "ACTUAL: " << actualResult << std::endl;
                //std::cout << "EXPECTED: " << expectedResult << std::endl;
            }
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

TEST_CASE("Cohesion Force between two Fluid Particles", "[Physics]") {
    Scenario scenario;
    scenario.COHESION = 0.25;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
            
            const Eigen::Vector3d actualResult = learnSPH::getCohesionForce(scenario, physicalData, fluidParticleIndex, neighborIndex);
            const Eigen::Vector3d expectedResult = fcohesion(scenario.COHESION, physicalData.fluidParticleMass, physicalData.fluidPositions[fluidParticleIndex], physicalData.fluidPositions[neighborIndex], 2.0*physicalData.h);

            if (actualResult.norm() > 0 || expectedResult.norm() > 0) {
                //std::cout << "ACTUAL: " << actualResult << std::endl;
                //std::cout << "EXPECTED: " << expectedResult << std::endl;
            }
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

TEST_CASE("Curvature Force between two Fluid Particles", "[Physics]") {
    Scenario scenario;
    scenario.COHESION = 0.25;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);
    learnSPH::ConfigData configData;
    configData.numberOfThreads = std::thread::hardware_concurrency();

    learnSPH::computeFluidSurfaceNormals(configData, physicalData, neighborhood);
    
    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
            
            const Eigen::Vector3d actualResult = learnSPH::getCurvatureForce(scenario, physicalData, fluidParticleIndex, neighborIndex);
            const Eigen::Vector3d expectedResult = fcurvature(scenario.COHESION, physicalData.fluidParticleMass, physicalData.fluidSurfaceNormals[fluidParticleIndex], physicalData.fluidSurfaceNormals[neighborIndex]); 

            if (actualResult.norm() > 0 || expectedResult.norm() > 0) {
                //std::cout << "ACTUAL: " << actualResult << std::endl;
                //std::cout << "EXPECTED: " << expectedResult << std::endl;
            }
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

TEST_CASE("Fluid Surface Normal Computation", "[Physics]") {
    Scenario scenario;
    scenario.COHESION = 0.25;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);
    learnSPH::ConfigData configData;
    configData.numberOfThreads = std::thread::hardware_concurrency();

    learnSPH::computeFluidSurfaceNormals(configData, physicalData, neighborhood);
    
    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        Eigen::Vector3d expectedNormal = Eigen::Vector3d(0, 0, 0);

        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
            expectedNormal += ni(2*physicalData.h, physicalData.fluidParticleMass, physicalData.densities[neighborIndex], physicalData.fluidPositions[fluidParticleIndex], physicalData.fluidPositions[neighborIndex], physicalData.h);
        });
        
        if (expectedNormal.norm() > 0 || physicalData.fluidSurfaceNormals[fluidParticleIndex].norm() > 0) {
            //std::cout << "ACTUAL: " << fluidSurfaceNormals[fluidParticleIndex] << std::endl;
            //std::cout << "EXPECTED: " << expectedNormal << std::endl;
        }
        REQUIRE(physicalData.fluidSurfaceNormals[fluidParticleIndex].x() == Approx(expectedNormal.x()));
        REQUIRE(physicalData.fluidSurfaceNormals[fluidParticleIndex].y() == Approx(expectedNormal.y()));
        REQUIRE(physicalData.fluidSurfaceNormals[fluidParticleIndex].z() == Approx(expectedNormal.z()));
    }
}

TEST_CASE("Final Cohesion Factor Computation", "[Physics]") {
    Scenario scenario;
    scenario.COHESION = 0.25;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.fluidPointSetId, [&](unsigned int neighborIndex) -> void {
            
            const double actualResult = learnSPH::getFinalCohesionFactor(scenario, physicalData, fluidParticleIndex, neighborIndex);
            const double expectedResult = (2 * scenario.FLUID_REST_DENSITY) / (physicalData.densities[fluidParticleIndex] + physicalData.densities[neighborIndex]);

            if (actualResult > 0 || expectedResult > 0) {
                //std::cout << "ACTUAL: " << actualResult << std::endl;
                //std::cout << "EXPECTED: " << expectedResult << std::endl;
            }
            REQUIRE(actualResult == Approx(expectedResult));
        });
    }
}

TEST_CASE("Adhesion Force between Fluid and Boundary Particle", "[Physics]") {
    Scenario scenario;
    scenario.COHESION = scenario.ADHESION = 0.25;
    learnSPH::PhysicalData physicalData = generateRandomPhysicalData(10000);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, 1000, 1);

    for (unsigned int fluidParticleIndex = 0; fluidParticleIndex < physicalData.fluidPositions.size(); fluidParticleIndex++) {
        neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, fluidParticleIndex, neighborhood.boundaryPointSetId, [&](unsigned int neighborIndex) -> void {
            
            const Eigen::Vector3d actualResult = learnSPH::getAdhesionForce(scenario, physicalData, fluidParticleIndex, neighborIndex);
            const Eigen::Vector3d expectedResult = fadhesion(scenario.ADHESION, physicalData.fluidParticleMass, physicalData.boundaryMasses[neighborIndex], physicalData.fluidPositions[fluidParticleIndex], physicalData.boundaryPositions[neighborIndex], 2*physicalData.h);

            if (actualResult.norm() > 0 || expectedResult.norm() > 0) {
                //std::cout << "ACTUAL: " << actualResult << std::endl;
                //std::cout << "EXPECTED: " << expectedResult << std::endl;
            }
            REQUIRE(actualResult.x() == Approx(expectedResult.x()));
            REQUIRE(actualResult.y() == Approx(expectedResult.y()));
            REQUIRE(actualResult.z() == Approx(expectedResult.z()));
        });
    }
}

// TODO: PBF Physics Tests
