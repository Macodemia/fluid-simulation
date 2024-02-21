#include "kernel.h"

static double kernelCubicFunctionLookupTable[learnSPH::kernel::LOOKUP_TABLE_SIZE];
static bool isLookupTableInitialized = false;
static double compactSupport = 0;

void learnSPH::kernel::computeLookupTables(const double h) {
    //if (!isLookupTableInitialized) {
        //isLookupTableInitialized = true;
        compactSupport = 2.0*h;
        for (int i = 0; i < LOOKUP_TABLE_SIZE; i++) {
            kernelCubicFunctionLookupTable[i] = learnSPH::kernel::kernelCubicFunction_computation(h, i * (compactSupport / (LOOKUP_TABLE_SIZE - 1)));
        }
    //}
}

// The q is not passed with division by h, instead for precomputation, h division is done inside the function
// Concrete: q = ||x1 - x2|| instead of ||x1 - x2|| / h
double learnSPH::kernel::cubicFunction(double q, const double h) {
    assert(q >= 0.0);
    assert(h > 0.0);
    q /= h;
    constexpr double alpha = 3.0 / (2.0 * M_PI);

    if (q < 1.0) {
        return alpha * (2.0 / 3.0 - pow(q, 2) + 0.5 * pow(q, 3));
    } else if (q < 2.0) {
        return alpha * (1.0 / 6.0 * pow(2.0 - q, 3));
    } else {
        return 0.0;
    }
}

double learnSPH::kernel::kernelCubicFunction_computation(const double h, const double q) {
    return (1.0 / pow(h, 3)) * cubicFunction(q, h);
}

double learnSPH::kernel::kernelCubicFunction(const double q) {
    const double lookupTableResolution = compactSupport / (LOOKUP_TABLE_SIZE - 1);
    const int index = (int)(q / lookupTableResolution);
    if (index >= LOOKUP_TABLE_SIZE) {
        return 0;
    }
    return kernelCubicFunctionLookupTable[index];
}

double learnSPH::kernel::kernelCubicFunction_computation(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double h) {
    const double q = (x1 - x2).norm();
    return learnSPH::kernel::kernelCubicFunction_computation(h, q);
}

double learnSPH::kernel::kernelCubicFunction(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) {
    const double q = (x1 - x2).norm();
    return kernelCubicFunction(q);
}

// The q is not passed with division by h, instead for precomputation, h division is done inside the function
// Concrete: q = ||x1 - x2|| instead of ||x1 - x2|| / h
double learnSPH::kernel::cubicGradFunction(double q, const double h) {
    assert(q >= 0.0);
    assert(h > 0.0);
    q /= h;
    constexpr double alpha = 3.0 / (2.0 * M_PI);

    if (q < 1.0) {
        return alpha * (-2.0 * q + 3.0 / 2.0 * (q * q));
    } else if (q < 2.0) {
        return alpha * (-1.0 / 2.0 * (2.0 - q) * (2.0 - q));
    } else {
        return 0.0;
    }
}

Eigen::Vector3d learnSPH::kernel::kernelGradCubicFunction(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double h) {
    if(x1 == x2){
        return Eigen::Vector3d (0, 0, 0);
    }

    const double q = (x1 - x2).norm();
    const double partWithoutVector = 1.0 / (h * h * h * h) * cubicGradFunction(q, h);
    const double partWithoutComponentSubtraction = partWithoutVector / (x1 - x2).norm();

    return Eigen::Vector3d(partWithoutComponentSubtraction * (x1.x() - x2.x()),
                           partWithoutComponentSubtraction * (x1.y() - x2.y()),
                           partWithoutComponentSubtraction * (x1.z() - x2.z()));
}

// Surface Tension Kernels

double learnSPH::kernel::cohesionKernel(const double r, const double c) {
    assert(r >= 0.0);
    assert(c >= 0.0);

    // Note: NOT using std::pow for performance reasons (when doing integer pows)
    const double cToPowerOfNine = c * c * c * c * c * c * c * c * c;
    const double cToPowerOfSix = c * c * c * c * c * c;

    const double rToPowerOfThree = r * r * r;
    const double cMinusRToPowerOfThree = (c - r) * (c - r) * (c - r);

    double returnValue;
    if (r <= c / 2.0) {
        returnValue = (32.0 / (M_PI * cToPowerOfNine)) * ((2 * cMinusRToPowerOfThree * rToPowerOfThree) - (cToPowerOfSix / 64.0));
    } else if (r <= c) {
        returnValue = (32.0 / (M_PI * cToPowerOfNine)) * (cMinusRToPowerOfThree * rToPowerOfThree);
    } else {
        returnValue = 0.0;
    }

    return returnValue;
}

double learnSPH::kernel::cohesionKernel(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double c) {
    const double r = (x1 - x2).norm();
    return cohesionKernel(r, c);
}

double learnSPH::kernel::adhesionKernel(const double r, const double c) {
    assert(r >= 0.0);
    assert(c >= 0.0);

    if (r >= c / 2.0 && r <= c) {
        return (0.007 / std::pow(c, 3.25)) * std::pow((-4 * r*r / c) + 6*r - 2*c, 1.0/4.0);
    } else {
        return 0;
    }
}

double learnSPH::kernel::adhesionKernel(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const double c) {
    const double r = (x1 - x2).norm();
    return adhesionKernel(r, c);
}
