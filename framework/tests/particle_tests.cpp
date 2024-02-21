#include "catch.hpp"

#include <thread>
#include <fstream>
#include <Eigen/Dense>
#include <limits>
#include "numeric"
#include "../learnSPH/kernel.h"
#include "../learnSPH/particle.h"
#include "../learnSPH/neighborhood.h"
#include "../learnSPH/data.h"

// NOTE: When simulating, this value is provided by the scenario. This is just defined here to avoid creating dummy scenarios.
static constexpr double FLUID_REST_DENSITY = 1000;
static constexpr double ETA = 1.2;

static learnSPH::PhysicalData getPhysicalData(const double width, const double height, const double depth, const int numberOfParticles, const Eigen::Vector3d& offset = {0, 0, 0}) {
    learnSPH::PhysicalData physicalData;
    learnSPH::particle::sampleParticles(physicalData, width, height, depth, FLUID_REST_DENSITY, numberOfParticles, offset);
    physicalData.computeH(FLUID_REST_DENSITY, ETA);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, FLUID_REST_DENSITY, std::thread::hardware_concurrency());
    return physicalData;
}

static learnSPH::PhysicalData getPhysicalData(const std::vector<Eigen::Vector3d>& boxPositions, const int numberOfParticles, const Eigen::Vector3d& offset = {0, 0, 0}) {
    learnSPH::PhysicalData physicalData;
    learnSPH::particle::sampleParticles(physicalData, boxPositions, FLUID_REST_DENSITY, numberOfParticles, offset);
    physicalData.computeH(FLUID_REST_DENSITY, ETA);
    learnSPH::Neighborhood neighborhood(physicalData.h, physicalData.fluidPositions, physicalData.boundaryPositions);
    physicalData.updateParticleDensitiesAndBoundaryVolumes(neighborhood, FLUID_REST_DENSITY, std::thread::hardware_concurrency());
    return physicalData;
}

TEST_CASE("Build Neighborhood - Particle should have itself as a neighbor", "[Particle]") {
    constexpr double h = 0.2;

    std::vector<Eigen::Vector3d> boundaryParticles;
    std::vector<Eigen::Vector3d> fluidParticles = {{1, 2, 3}};

    learnSPH::Neighborhood neighborhood(h, fluidParticles, boundaryParticles);

    bool hadItselfAsNeighbor = false;
    neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, 0, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
        if (neighborIndex == 0) {
            hadItselfAsNeighbor = true;
        }
        REQUIRE(neighborIndex == 0); // To make sure only one neighbor (itself) was found
    });
    REQUIRE(hadItselfAsNeighbor == true);

    const std::vector<double> hValues = { 0.0000000001, 0.2, 1000.0 };

    for (const double hValue : hValues) {
        // std::get<0> because we only need the particles and not the massPerParticle
        fluidParticles = std::get<0>(learnSPH::particle::sampleParticles(learnSPH::particle::UNIT_CUBE, FLUID_REST_DENSITY, 1000));
        neighborhood = learnSPH::Neighborhood(hValue, fluidParticles, boundaryParticles);

        for (int particleIndex = 0; particleIndex < fluidParticles.size(); particleIndex++) {
            bool foundItselfAsNeighbor = false;
            neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, particleIndex, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
                if (neighborIndex == particleIndex) {
                    foundItselfAsNeighbor = true;
                }
            });
            REQUIRE(foundItselfAsNeighbor == true);
        }
    }
}

// TODO: !!!
#if 0
TEST_CASE("Build Neighborhood - Should compute correct neighbors", "[Particle]") {
    constexpr double h = 0.25;
    const std::vector<std::pair<std::vector<Eigen::Vector3d>, std::vector<std::vector<int>>>> inputExpectedResultPairs = {
            {
                // Input Particles - 3 Clusters that are completely separate
                {
                    {0, 0, 0}, {0.1, 0.1, 0.1}, {0, 0.2, 0},
                    {2, 2, 2}, {2.1, 1.9, 2},
                    {5, 5, 5 }
                },
                // Expected Results
                {
                    { 0, 1, 2 }, { 1, 0, 2 }, { 2, 0, 1},
                    { 3, 4 }, { 4, 3 },
                    { 5 }
                }
            },
            {
                // Input Particles - One middle particle is in the neighborhood of two edge particles, the edge particles are not in the neighborhood of each other
                {
                    { 0, 0, 0 }, { 0.25, 0.25, 0.25 }, { 0.4, 0.4, 0.4 }
                },
                // Expected Results
                {
                    { 0, 1 }, { 1, 0, 2 }, { 2, 1 }
                }
            },
            {
                // Input Particles - All particles are neighbours of each other
                {
                    { 0, 0, 0 }, { 0.1, 0.1, 0.1 }, { 0.2, 0.2, 0.2 }, { 0.25, 0.25, 0.25 }
                },
                // Expected Results
                {
                    { 0, 1, 2, 3 }, { 1, 0, 2, 3 }, { 2, 0, 1, 3 }, { 3, 0, 1, 2 }
                }
            }
    };

    for (auto const& inputExpectedResultPair : inputExpectedResultPairs) {
        std::vector<Eigen::Vector3d> particles = inputExpectedResultPair.first;
        std::vector<std::vector<int>> expectedNeighbors = inputExpectedResultPair.second;

        std::vector<std::vector<int>> resultingNeighbors = learnSPH::particle::buildNeighborhood(particles, h);

        REQUIRE(expectedNeighbors == resultingNeighbors);
    }
}
#endif

TEST_CASE("Sample Particles - All particles should be inside the box", "[Particle]") {
    // inputs: width, height, depth
    std::vector<std::tuple<double, double, double>> inputs = {
            { 1, 1, 1 }, { 0.1, 0.2, 0.3 }, { 10, 0.1, 0.1 }, { sqrt(2), sqrt(3), sqrt(M_PI) }
    };

    for (std::tuple<double, double, double> input : inputs) {
        double width, height, depth;
        std::tie(width, height, depth) = input;

        std::vector<Eigen::Vector3d> particles;
        double mass;
        double samplingDistance;
        std::tie(particles, mass, samplingDistance) = learnSPH::particle::sampleParticles(width, height, depth, FLUID_REST_DENSITY, 1000);

        for (const Eigen::Vector3d& particle : particles) {
            REQUIRE(particle.x() >= 0);
            REQUIRE(particle.x() <= Approx(width));
            REQUIRE(particle.y() >= 0);
            REQUIRE(particle.y() <= Approx(height));
            REQUIRE(particle.z() >= 0);
            REQUIRE(particle.z() <= Approx(depth));
        }
    }
}

TEST_CASE("Sample Particles - Box should be filled with particles", "[Particle]") {
    // inputs: width, height, depth
    std::vector<std::tuple<double, double, double>> inputs = {
            { 1, 1, 1 }, { 0.1, 0.2, 0.3 }, { 10, 0.1, 0.1 }, { sqrt(2), sqrt(3), sqrt(M_PI) }
    };

    for (std::tuple<double, double, double> input : inputs) {
        double width, height, depth;
        std::tie(width, height, depth) = input;

        std::vector<Eigen::Vector3d> particles;
        double mass;
        double samplingDistance;
        std::tie(particles, mass, samplingDistance) = learnSPH::particle::sampleParticles(width, height, depth, FLUID_REST_DENSITY, 10000);

        // Calculate smallest bounding box of sampled particles
        double minX = std::numeric_limits<double>::max(), maxX = -std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max(), maxY = -std::numeric_limits<double>::max();
        double minZ = std::numeric_limits<double>::max(), maxZ = -std::numeric_limits<double>::max();

        for (const Eigen::Vector3d& particle : particles) {
            minX = std::min(particle.x(), minX);
            maxX = std::max(particle.x(), maxX);
            minY = std::min(particle.y(), minY);
            maxY = std::max(particle.y(), maxY);
            minZ = std::min(particle.z(), minZ);
            maxZ = std::max(particle.z(), maxZ);
        }

        Approx targetZero = Approx(0).margin(0.1);
        REQUIRE(minX == targetZero);
        REQUIRE(minY == targetZero);
        REQUIRE(minZ == targetZero);

        Approx targetWidth = Approx(width).margin(0.1);
        Approx targetHeight = Approx(height).margin(0.1);
        Approx targetDepth = Approx(depth).margin(0.1);
        REQUIRE(maxX == targetWidth);
        REQUIRE(maxY == targetHeight);
        REQUIRE(maxZ == targetDepth);
    }
}

TEST_CASE("Sample Particles - All particles in a neighborhood should have the same average distance", "[Particle]") {
    constexpr double h = 0.05;
    // inputs: width, height, depth
    std::vector<std::tuple<double, double, double>> inputs = {
            { 1, 1, 1 }, { 0.1, 0.2, 0.3 }, { 10, 0.1, 0.1 }, { sqrt(2), sqrt(3), sqrt(M_PI) }
    };

    for (std::tuple<double, double, double> input : inputs) {
        double width, height, depth;
        std::tie(width, height, depth) = input;

        std::vector<Eigen::Vector3d> particles;
        std::vector<Eigen::Vector3d> boundaryParticles;
        double mass;
        double samplingDistance;
        std::tie(particles, mass, samplingDistance) = learnSPH::particle::sampleParticles(width, height, depth, FLUID_REST_DENSITY, 5000);

        learnSPH::Neighborhood neighborhood(h, particles, boundaryParticles);

        std::vector<double> meanDistances;
        meanDistances.reserve(particles.size());

        for (int i = 0; i < particles.size(); i++) {
            double meanDistance = 0;
            int numberOfNeighbors = 0;
            neighborhood.forEachNeighborDo(neighborhood.fluidPointSetId, i, neighborhood.fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
                meanDistance += (particles[i] - particles[neighborIndex]).norm();
                numberOfNeighbors++;
            });

            meanDistance /= numberOfNeighbors;

            meanDistances.push_back(meanDistance);
        }

        for (int i = 0; i < meanDistances.size(); i++) {
            for (int j = i; j < meanDistances.size(); j++) {
                Approx targetDistance = Approx(meanDistances[j]).epsilon(0.5);
                REQUIRE(meanDistances[i] == targetDistance);
            }
        }
    }
}

TEST_CASE("Density Estimation - Empty Point Set", "[Particle]") {
    learnSPH::PhysicalData physicalData = getPhysicalData(learnSPH::particle::UNIT_CUBE, 0);

    REQUIRE(physicalData.densities.empty());
}

TEST_CASE("Density Estimation - Single Point Set", "[Particle]") {
    learnSPH::PhysicalData physicalData = getPhysicalData(learnSPH::particle::UNIT_CUBE, 1);

    REQUIRE(physicalData.densities.size() == 1);
    REQUIRE(physicalData.densities[0] == Approx(183.654379).epsilon(0.01));

    // Manual Result Computation (For unit cube and one particle):
    // Volume of Cuboid = 1^3 = 1
    // Fluid Rest Density = 997 kg/m^3
    // Mass Per Particle = 997 kg/m^3 * 1 m^3 = 997 kg
    // h = ETA * representativeDiameter = 1.2 * pow(Mass Per Particle / Fluid Rest Density, 1/3) = 1.2 * pow(997 / 997, 1/3) = 1.2 * 1 = 1.2
    // Kernel Weight for Single Particle = (1 / h^3) * (1 / PI) = (1 / 1.2^3) * (1 / PI) = 0.184207
    // Density = Mass Per Particle * Kernel Weight = 997 * 0.184207 = 183.654379
}

TEST_CASE("Density Estimation - Average Point Density Values", "[Particle]") {
    learnSPH::PhysicalData physicalData = getPhysicalData(learnSPH::particle::UNIT_CUBE, 10000);

    /*
    std::ofstream out;
    out.open("points.txt");
    for (const Eigen::Vector3d& particle : fluidParticles) {
        out << particle.x() << " " << particle.y() << " " << particle.z() << std::endl;
    }
    out.close();

    std::ofstream out2;
    out2.open("points_not_rest.txt");*/
    
    Approx targetDensity = Approx(FLUID_REST_DENSITY).epsilon(0.01);
    
    for (int i = 0; i < physicalData.fluidPositions.size(); i++) {
        const Eigen::Vector3d& particle = physicalData.fluidPositions[i];

        if (physicalData.densities[i] != targetDensity) {
            //out2 << particle.x() << " " << particle.y() << " " << particle.z() << std::endl;
           
            // 0 and 1 because we are using learnSPH::particle::UNIT_CUBE 
            const double minDistanceToXBoundary = std::min(particle.x() - 0, 1 - particle.x());
            const double minDistanceToYBoundary = std::min(particle.y() - 0, 1 - particle.y());
            const double minDistanceToZBoundary = std::min(particle.z() - 0, 1 - particle.z());
            const double minDistanceToAnyBoundary = std::min(minDistanceToXBoundary, std::min(minDistanceToYBoundary, minDistanceToZBoundary));
           
            //std::cout << "Min Distance To Any Boundary: " << minDistanceToAnyBoundary << " < " << h << std::endl;

            REQUIRE(minDistanceToAnyBoundary < 2 * physicalData.h); // 2h <=> compact support
        }
    }
    //out2.close();
}
