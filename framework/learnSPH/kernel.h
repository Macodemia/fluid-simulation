#pragma once

#include <cassert>
#include <cmath>
#include <Eigen/Dense>

#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    namespace kernel {
        // TODO: Revise Kernel Precomputation
        // SPH Kernel
        constexpr int LOOKUP_TABLE_SIZE = 1000000;

        void computeLookupTables(const double h);

        // The q is not passed with division by h, instead for precomputation, h division is done inside the function
        // Concrete: q = ||x1 - x2|| instead of ||x1 - x2|| / h
        double cubicFunction(double q, const double h);

        double kernelCubicFunction_computation(const double h, const double q);
        double kernelCubicFunction(const double q);

        double kernelCubicFunction_computation(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double h);
        double kernelCubicFunction(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2);

        // The q is not passed with division by h, instead for precomputation, h division is done inside the function
        // Concrete: q = ||x1 - x2|| instead of ||x1 - x2|| / h
        double cubicGradFunction(double q, const double h);

        Eigen::Vector3d kernelGradCubicFunction(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double h);

        // Surface Tension Kernels
        double cohesionKernel(const double r, const double c);
        double cohesionKernel(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double c);

        double adhesionKernel(const double r, const double c);
        double adhesionKernel(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double c);
    };
};
