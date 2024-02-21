#pragma once
#include <CompactNSearch/CompactNSearch.h>
#include <Eigen/Dense>

namespace learnSPH {
    struct Neighborhood {
        CompactNSearch::NeighborhoodSearch nsearch;
        unsigned int fluidPointSetId;
        unsigned int boundaryPointSetId;

        //Note: 2 * h because that is our compact support. If the kernel function is changed, this needs to be updated
        Neighborhood(const double h, const std::vector<Eigen::Vector3d>& fluidParticles, std::vector<Eigen::Vector3d>& boundaryParticles) : nsearch(2 * h) {
            // CAREFUL: CompactNSearch stores the given double pointer -> If it is invalid (e.g. vector push_back or resize), it needs to be recomputed!
            fluidPointSetId = nsearch.add_point_set((double*)&fluidParticles[0], fluidParticles.size());
            boundaryPointSetId = nsearch.add_point_set((double*)&boundaryParticles[0], boundaryParticles.size());
            update();
        }

        inline void update() {
            nsearch.find_neighbors();
        }

        template<typename Function>
        void forEachNeighborDo(const unsigned int pointSetId, const unsigned int particleIndex, const unsigned int pointSetIdOfNeighbors, Function callback) const {
            const CompactNSearch::PointSet& pointSet = nsearch.point_set(pointSetId);
            // Should have itself as a neighbor (if both pointSetIds are the same, because when a fluid particle looks at boundary particle neighbors, it cannot see itself)
            if (pointSetId == pointSetIdOfNeighbors) {
                callback(particleIndex);
            }

            // Remaining neighbors (without itself)
            for (size_t neighborNSearchIndex = 0; neighborNSearchIndex < pointSet.n_neighbors(pointSetIdOfNeighbors, particleIndex); neighborNSearchIndex++) {
                const unsigned int neighborIndex = pointSet.neighbor(pointSetIdOfNeighbors, particleIndex, neighborNSearchIndex);
                callback(neighborIndex);
            }
        }
    };
}
