#include "catch.hpp"

#include <cmath>
#include <iostream>
#include <tuple>
#include <random>
#include <limits>
#include <array>

#include "../learnSPH/kernel.h"

static double generateRandomDouble(double min, double max) {
    std::uniform_real_distribution<double> distribution(min, max);
    std::random_device r;
    std::default_random_engine randomEngine(r());
    return distribution(randomEngine);
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

TEST_CASE("Cubic Spline", "[Kernel]") {
    // first entry is param, second entry is expected result
    const std::tuple<double, double> values[] = {
            {0,    1.0 / M_PI},
            {0.5,  0.228785},
            {1.0,  1.0 / (4 * M_PI)},
            {1.5,  0.00994718},
            {2.0,  0},
            {2.5,  0},
            {1000, 0}
    };

    for (auto value : values) {
        const double param = std::get<0>(value);
        const double expectedResult = std::get<1>(value);

        const double actualResult = learnSPH::kernel::cubicFunction(param, 1.0); // TODO: Maybe multiple h values?

        REQUIRE(actualResult == Approx(expectedResult));
    }
}

TEST_CASE("Cubic Spline compared with generated function", "[Kernel]") {
    for (int i = 0; i < 10000; i++) {
        const double q = generateRandomDouble(0.0, 10.0);
        const double expectedResult = s(q);
        const double actualResult = learnSPH::kernel::cubicFunction(q, 1.0);

        REQUIRE(actualResult == Approx(expectedResult));
    }    
}

TEST_CASE("Kernel (Cubic) Function", "[Kernel]") {
    // first entry is vec1, second entry is vec2, third entry is h, fourth entry is expected result
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, double, double> values[] = {
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 1), 0.1,  318.3098861837906},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 1), 0.5,  2.546479089470325},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 1), 1.5,  0.09431404035075278},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 0.1,  0},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 0.5,  0},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 1.0,  0.07957747154594767},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 1.5,  0.05239668908375155},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 10.0, 0.00031377397030567163},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(5, 5, 5), 0.5,  0},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(5, 5, 5), 8,    0.0002251416329382914},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(5, 5, 5), 3.5,  1.60212569114751e-08},
    };

    for (auto value : values) {
        const Eigen::Vector3d vec1 = std::get<0>(value);
        const Eigen::Vector3d vec2 = std::get<1>(value);
        const double h = std::get<2>(value);
        const double expectedResult = std::get<3>(value);

        const double actualResult = learnSPH::kernel::kernelCubicFunction_computation(vec1, vec2, h);

        REQUIRE(actualResult == Approx(expectedResult));
    }
}

TEST_CASE("Kernel (Cubic) Function Precomputation", "[Kernel]") {
    // first entry is vec1, second entry is vec2, third entry is h
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, double> values[] = {
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 1), 0.1},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 1), 0.5},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 1), 1.5},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 0.1},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 0.5},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 1.0},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 1.5},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(1, 1, 0), 10.},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(5, 5, 5), 0.5},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(5, 5, 5), 8},
            {Eigen::Vector3d(1, 1, 1), Eigen::Vector3d(5, 5, 5), 3.5},
    };

    for (auto value : values) {
        const Eigen::Vector3d vec1 = std::get<0>(value);
        const Eigen::Vector3d vec2 = std::get<1>(value);
        const double h = std::get<2>(value);
    
        learnSPH::kernel::computeLookupTables(h);

        const double recomputedResult = learnSPH::kernel::kernelCubicFunction_computation(vec1, vec2, h);
        const double precomputedResult = learnSPH::kernel::kernelCubicFunction(vec1, vec2);

        std::cout << "Expected: " << recomputedResult << " - Actual: " << precomputedResult << std::endl;

        if (precomputedResult == 0.0 || recomputedResult == 0.0) {
            REQUIRE(recomputedResult == Approx(precomputedResult).margin(0.000000000000000000001));
        } else {
            REQUIRE(recomputedResult == Approx(precomputedResult).epsilon(0.001));
        }
    }
}

// TODO: For each test case: write down what could go wrong with random generated 3d vectors and what we do to mitigate that
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

// TODO: How to do the Unity property (Integral)

static double hValues[] = {0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 10.0};

TEST_CASE("Kernel (Cubic) Function should have Symmetry property", "[Kernel]") {
    std::vector<Eigen::Vector3d> vecs = generateRandomVectors3d(100);
    for (double h : hValues) {
        for (const Eigen::Vector3d& vec1 : vecs) {
            for (const Eigen::Vector3d& vec2 : vecs) {
                const double result1 = learnSPH::kernel::kernelCubicFunction_computation(vec1, vec2, h);
                const double result2 = learnSPH::kernel::kernelCubicFunction_computation(vec2, vec1, h);
                REQUIRE(result1 == Approx(result2));
            }
        }
    }
}

TEST_CASE("Kernel (Cubic) Function should have Delta property", "[Kernel]") {
    std::vector<Eigen::Vector3d> vecs = generateRandomVectors3d(100);
    constexpr double h = 1e-8; // because the property uses lim h -> 0

    for (const Eigen::Vector3d& vec1 : vecs) {
        for (const Eigen::Vector3d& vec2 : vecs) {
            const double result1 = learnSPH::kernel::kernelCubicFunction_computation(vec1, vec2, h);

            if (vec1 == vec2) {
                REQUIRE(result1 == Approx(std::numeric_limits<double>::infinity()));
            } else {
                REQUIRE(result1 == Approx(0));
            }
        }
    }
}

TEST_CASE("Kernel (Cubic) Function should have Non-negative property", "[Kernel]") {
    std::vector<Eigen::Vector3d> vecs = generateRandomVectors3d(100);

    for (double h : hValues) {
        for (const Eigen::Vector3d& vec1 : vecs) {
            for (const Eigen::Vector3d& vec2 : vecs) {
                const double result1 = learnSPH::kernel::kernelCubicFunction_computation(vec1, vec2, h);
                REQUIRE(result1 >= 0);
            }
        }
    }
}

TEST_CASE("Kernel (Cubic) Function should have Compactness property", "[Kernel]") {
    std::vector<Eigen::Vector3d> vecs = generateRandomVectors3d(100);

    // Add two vectors that are very close together, so that at least one vector pair hits the else case below
    vecs.emplace_back(1.0, 1.0, 1.0);
    vecs.emplace_back(1.01, 1.0, 1.0);

    for (double h : hValues) {
        for (const Eigen::Vector3d& vec1 : vecs) {
            for (const Eigen::Vector3d& vec2 : vecs) {
                const double result1 = learnSPH::kernel::kernelCubicFunction_computation(vec1, vec2, h);
                const double norm = (vec1 - vec2).stableNorm();

                if (norm > 2.0 * h) {
                    REQUIRE(result1 == Approx(0));
                }
            }
        }
    }
}

static Eigen::Vector3d finiteDifferenceApproximation(const double epsilon, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2,
                                                     const double h) {
    using namespace learnSPH::kernel;

    const double qPlusX = (x1 - x2 + epsilon * Eigen::Vector3d(1, 0, 0)).norm();
    const double qMinusX = (x1 - x2 - epsilon * Eigen::Vector3d(1, 0, 0)).norm();

    const double qPlusY = (x1 - x2 + epsilon * Eigen::Vector3d(0, 1, 0)).norm();
    const double qMinusY = (x1 - x2 - epsilon * Eigen::Vector3d(0, 1, 0)).norm();

    const double qPlusZ = (x1 - x2 + epsilon * Eigen::Vector3d(0, 0, 1)).norm();
    const double qMinusZ = (x1 - x2 - epsilon * Eigen::Vector3d(0, 0, 1)).norm();

    return Eigen::Vector3d((kernelCubicFunction_computation(h, qPlusX) - kernelCubicFunction_computation(h, qMinusX)) / (2.0 * epsilon),
                           (kernelCubicFunction_computation(h, qPlusY) - kernelCubicFunction_computation(h, qMinusY)) / (2.0 * epsilon),
                           (kernelCubicFunction_computation(h, qPlusZ) - kernelCubicFunction_computation(h, qMinusZ)) / (2.0 * epsilon));
}

static double relativeDifference(double a, double b) {
    // https://en.wikipedia.org/wiki/Relative_change_and_difference
    // (Relative Difference with max(|x|, |y|) for function f(x,y)
    // NOTE: WHEN BOTH APPROX AND ANALYTICAL ARE 0 => DIVISION BY ZERO
    // Smaller is better
    // NOTE: Can get over 100% (of course: Difference can be larger than "percentage basis")

    if (a != 0 || b != 0) {
        return 100 * abs(a - b) / (fmax(abs(a), abs(b)));
    }
    return 0;
}

TEST_CASE("Analytical gradient should be close to finite difference approximation", "[Kernel]") {
    const double epsilon = 1e-6;

    for (double h : hValues) {
        std::vector<Eigen::Vector3d> vecs = generateRandomVectors3d(100);

        for (const Eigen::Vector3d& vec1 : vecs) {
            for (const Eigen::Vector3d& vec2 : vecs) {
                if (vec1 == vec2) { continue; }

                const Eigen::Vector3d approx = finiteDifferenceApproximation(epsilon, vec1, vec2, h);
                const Eigen::Vector3d analytical = learnSPH::kernel::kernelGradCubicFunction(vec1, vec2, h);

                // NOTE
                // We do not want to look at absolute differences as the result will look good if both values are very small
                // e.g. if approx.x is 1.86024e-06 and analytical.x is 7.93941e-06 then the absolute difference is "only" 6.07917e-06
                // however the values should still not be considered similar
                // Therefore: Use a relative difference calculation
                const Eigen::Vector3d relativeDifferenceVector = Eigen::Vector3d(relativeDifference(approx.x(), analytical.x()),
                                                                                 relativeDifference(approx.y(), analytical.y()),
                                                                                 relativeDifference(approx.z(), analytical.z()));

                const double maxRelativeDifference = fmax(fmax(relativeDifferenceVector.x(), relativeDifferenceVector.y()), relativeDifferenceVector.z());

                /*std::cout << "\nAnalytical gradient should be close to finite difference approximation" << std::endl;
                std::cout << "Approx:\n" << approx << std::endl << "Analytical:\n" << analytical << std::endl
                          << "MaxRelativeDiff:" << maxRelativeDifference
                          << std::endl;*/

                REQUIRE(maxRelativeDifference < 1.0); // 1.0 -> 1% Relative Difference
            }
        }
    }
}

TEST_CASE("Cohesion Kernel should compute correct values", "[Kernel]") {
    const double h = 1.0; // c = 2.0

    // first entry is r, second entry is expected result
    const std::tuple<double, double> values[] = {
            {0, -0.019894367},
            {0.5, -0.00310849},
            {0.9, 0.0187126},
            {1.0, 0.0198944},
            {1.1, 0.0193035},
            {1.5, 0.00839294},
            {1.9, 0.000136455},
            {2.0, 0},
            {2.1, 0},
            {3.0, 0},
            {50.0, 0}
    };

    for (auto value : values) {
        const double r = std::get<0>(value);
        const double expectedResult = std::get<1>(value);

        const double actualResult = learnSPH::kernel::cohesionKernel(r, 2 * h);

        REQUIRE(actualResult == Approx(expectedResult));
    }
}

TEST_CASE("Adhesion Kernel should compute correct values", "[Kernel]") {
    const double h = 1.0; // c = 2.0

    // first entry is r, second entry is expected result
    const std::tuple<double, double> values[] = {
            {0, 0},
            {0.5, 0},
            {0.9, 0},
            {1.0, 0},
            {1.1, 0.000479257},
            {1.5, 0.000618718},
            {1.9, 0.000479257},
            {2.0, 0},
            {2.1, 0},
            {2.5, 0},
            {3.0, 0},
            {50.0, 0}
    };

    for (auto value : values) {
        const double r = std::get<0>(value);
        const double expectedResult = std::get<1>(value);

        const double actualResult = learnSPH::kernel::adhesionKernel(r, 2 * h);

        REQUIRE(actualResult == Approx(expectedResult));
    }
}
