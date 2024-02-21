#pragma once

#include <Eigen/Dense>
#include <thread>
#include "data.h"

namespace learnSPH {
    namespace particle {
        const std::vector<Eigen::Vector3d> UNIT_CUBE = {{0, 1, 0},
                                                        {1, 1, 0},
                                                        {1, 1, 1},
                                                        {0, 1, 1},
                                                        {0, 0, 0},
                                                        {1, 0, 0},
                                                        {1, 0, 1},
                                                        {0, 0, 1}};

        // Return tuple is: (sampledParticles, calculatedMassPerParticle)
        std::tuple<std::vector<Eigen::Vector3d>, double, double> sampleParticles(const double width, const double height, const double depth, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset = {0, 0, 0}, const Eigen::Vector3d& rotation = {0, 0, 0});
        std::tuple<std::vector<Eigen::Vector3d>, double, double> sampleParticles(const std::vector<Eigen::Vector3d>& boxPositions, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset = {0, 0, 0}, const Eigen::Vector3d& rotation = {0, 0, 0});
        void sampleParticles(learnSPH::PhysicalData& physicalData, const std::vector<Eigen::Vector3d>& boxPositions, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset = {0, 0, 0}, const Eigen::Vector3d& rotation = {0, 0, 0});
        void sampleParticles(learnSPH::PhysicalData& physicalData, const double width, const double height, const double depth, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset = {0, 0, 0}, const Eigen::Vector3d& rotation = {0, 0, 0});

        std::vector<Eigen::Vector3d> sampleTriangle(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C, const double samplingDistance, const bool hexagonalSampling = true);

        std::vector<Eigen::Vector3d> sampleBoundaryBox(const double width, const double height, const double depth, const double samplingDistance, Eigen::Vector3d offset = {0, 0, 0});

        std::vector<Eigen::Vector3d> sampleBoundaryMesh(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::vector<int>>& faces, const double samplingDistance);

        std::vector<Eigen::Vector3d> sampleObjFile(const std::string& file, const double samplingDistance, const Eigen::Vector3d& offset = Eigen::Vector3d(0, 0, 0), const double scale = 1.0);

    };
};
