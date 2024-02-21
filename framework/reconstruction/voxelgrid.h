#pragma once
#include <thread>
#include <memory>
#include "../learnSPH/util/marching_cubes_lut.hpp"
#include "../learnSPH/kernel.h"
#include "../learnSPH/util/vtk_writer.h"
#include "../learnSPH/particle.h"
#include "../learnSPH/neighborhood.h"

typedef std::array<uint64_t, 8> GridCell; // 8 Vertex indices

enum EdgeIndexDirection {
    X = 0,
    Y = 1,
    Z = 2
};

struct ProgressData {
    int particleCount;
    int particleIndex;
    std::string frameFile;
    int endIndex;
    bool threadFinished = false;
};

struct Voxelgrid {
    std::string gridVerticesFileName;
    std::vector<double> levelSetValues;
    Eigen::Vector3d offset;
    double cellEdgeLength{};
    uint64_t cellsX{}, cellsY{}, cellsZ{};
    uint64_t verticesX{};
    uint64_t verticesY{};
    uint64_t verticesZ{};

    Voxelgrid() {
        levelSetValues.reserve(300000);
    }

    void init(const Eigen::Vector3d& offset, const double boundingBoxWidth, const double boundingBoxHeight, const double boundingBoxDepth, const double cellEdgeLength, const std::string gridVerticesFileName) {
        this->offset = offset;
        this->cellEdgeLength = cellEdgeLength;
        std::tie(this->cellsX, this->cellsY, this->cellsZ) = getCellCount(boundingBoxWidth, boundingBoxHeight, boundingBoxDepth, cellEdgeLength);
        this->verticesX = cellsX + 1;
        this->verticesY = cellsY + 1;
        this->verticesZ = cellsZ + 1;
        this->gridVerticesFileName = gridVerticesFileName;
    }

    inline uint64_t getVertexCount() const {
        return verticesX * verticesY * verticesZ;
    }

    inline Eigen::Vector3d getVertexCoordinates(const uint64_t i, const uint64_t j, const uint64_t k) const {
        return { offset.x() + i * cellEdgeLength, offset.y() + j * cellEdgeLength, offset.z() + k * cellEdgeLength};
    }

    [[nodiscard]] uint64_t getVertexIndex(const uint64_t i, const uint64_t j, const uint64_t k) const {
        return (i * verticesY * verticesZ) + (j * verticesZ) + k;
    }

    [[nodiscard]] uint64_t getCellIndex(const uint64_t i, const uint64_t j, const uint64_t k) const {
        return (i * cellsY * cellsZ) + (j * cellsZ) + k;
    }

    [[nodiscard]] static uint64_t getEdgeIndex(const uint64_t indexOfVertexOfOrigin, EdgeIndexDirection dir) {
        return (3 * indexOfVertexOfOrigin) + dir;
    }

    [[nodiscard]] std::tuple<uint64_t, uint64_t, uint64_t> getijkFromV(const uint64_t V) const {
        uint64_t k = V % verticesZ;
        uint64_t j = ((V - k) / verticesZ) % verticesY;
        uint64_t i = (((V - k) / verticesZ) - j) / verticesY;
        return std::make_tuple(i, j, k);
    }

    [[nodiscard]] bool isVAtBorderInDir(const uint64_t V, EdgeIndexDirection dir) const {
        uint64_t i, j, k;
        std::tie(i, j, k) = getijkFromV(V);

        switch (dir) {
            case EdgeIndexDirection::X:
                return i < verticesX - 1;
            case EdgeIndexDirection::Y:
                return j < verticesY - 1;
            case EdgeIndexDirection::Z:
                return k < verticesZ - 1;
        }
    }

    [[nodiscard]] std::tuple<std::vector<Eigen::Vector3d>, std::vector<std::array<int, 3>>> marchingCubes() const {
        // 3.
        // intersectionPoints contains 3d vectors that represent the points from which triangles should be created
        std::vector<Eigen::Vector3d> intersectionPoints;
        // edgeIntersectionPointMapping is a mapping of global edge identifiers to intersection point indices
        std::unordered_map<uint64_t, uint64_t> edgeIntersectionPointMapping;

        for (int V = 0; V < getVertexCount(); V++) {
            for (int dir = 0; dir < 3; dir++) {
                if (!isVAtBorderInDir(V, (EdgeIndexDirection)dir)) {
                    continue;
                }
                uint64_t i, j, k;
                std::tie(i, j, k) = getijkFromV(V);

                const uint64_t Va = getVertexIndex(i, j, k);
                uint64_t Vb;

                switch (dir) {
                    case EdgeIndexDirection::X:
                        Vb = getVertexIndex(i+1, j, k);
                        break;
                    case EdgeIndexDirection::Y:
                        Vb = getVertexIndex(i, j+1, k);
                        break;
                    case EdgeIndexDirection::Z:
                        Vb = getVertexIndex(i, j, k+1);
                        break;
                    default:
                        std::cout << "Error: Checking for non-existing dir (" << dir << "). Stopping marching cubes!" << std::endl;
                        std::vector<std::array<int, 3>> mesh;
                        return std::make_tuple(intersectionPoints, mesh);
                }

                const Eigen::Vector3d xa = getVertexCoordinates(i, j, k);
                uint64_t ib, jb, kb;
                std::tie(ib, jb, kb) = getijkFromV(Vb);
                const Eigen::Vector3d xb = getVertexCoordinates(ib, jb, kb);

                const double capitalPhi_xa = levelSetValues[Va];
                const double capitalPhi_xb = levelSetValues[Vb];

                if (capitalPhi_xa < 0 && capitalPhi_xb > 0 || capitalPhi_xa > 0 && capitalPhi_xb < 0) {
                    const double alpha = (capitalPhi_xa) / (-capitalPhi_xb + capitalPhi_xa);
                    const Eigen::Vector3d xs = (1.0 - alpha) * xa + alpha * xb;
                    intersectionPoints.push_back(xs);
                    const uint64_t E = getEdgeIndex(Va, (EdgeIndexDirection)dir);
                    edgeIntersectionPointMapping[E] = intersectionPoints.size() - 1;
                }
            }
        }

        constexpr uint64_t localEdgeToVertexAndDirMap[12][2] = {
                {0, EdgeIndexDirection::X}, // 0
                {1, EdgeIndexDirection::Y}, // 1
                {3, EdgeIndexDirection::X}, // 2
                {0, EdgeIndexDirection::Y}, // 3
                {4, EdgeIndexDirection::X}, // 4
                {5, EdgeIndexDirection::Y}, // 5
                {7, EdgeIndexDirection::X}, // 6
                {4, EdgeIndexDirection::Y}, // 7
                {0, EdgeIndexDirection::Z}, // 8
                {1, EdgeIndexDirection::Z}, // 9
                {2, EdgeIndexDirection::Z}, // 10
                {3, EdgeIndexDirection::Z}, // 11
        };

        // 4.
        std::vector<std::array<int, 3>> mesh;
        for (uint64_t i = 0; i < cellsX; i++) {
            for (uint64_t j = 0; j < cellsY; j++) {
                for (uint64_t k = 0; k < cellsZ; k++) {
                    GridCell gridCell = getGridCell(i, j, k);

                    std::array<bool, 8> vertexSigns{};
                    for (int j = 0; j < gridCell.size(); j++) {
                        vertexSigns[j] = levelSetValues[gridCell[j]] > 0;
                    }

                    std::array<std::array<int, 3>, 5> triangles = getMarchingCubesCellTriangulation(vertexSigns);
                    for (const std::array<int, 3>& localEdgeTriangle : triangles) {
                        if (localEdgeTriangle[0] == -1) { break; } // => We have already saved all needed triangles

                        // 3 intersection point indices
                        std::array<int, 3> meshTriangle{};

                        for (int i = 0; i < 3; i++) {
                            const uint64_t V = gridCell[localEdgeToVertexAndDirMap[localEdgeTriangle[i]][0]];
                            const auto dir = (EdgeIndexDirection)localEdgeToVertexAndDirMap[localEdgeTriangle[i]][1];
                            const uint64_t E = getEdgeIndex(V, dir);
                            const auto got = edgeIntersectionPointMapping.find(E);
                            assert(got != edgeIntersectionPointMapping.end());
                            meshTriangle[i] = got->second;
                        }

                        mesh.push_back(meshTriangle);
                    }
                }
            }
        }

        // mesh now contains a list of triangles (3 intersection point indices)
        return std::make_tuple(intersectionPoints, mesh);
    }

    void getLevelSetValues_SPH(const std::vector<Eigen::Vector3d>& fluidParticles, const learnSPH::Neighborhood& neighborhood, const double h, const double c, ProgressData* const progress, int maxAllowedThreadSplits = 0) {
        // This function only handles threading of the level set value computation.
        // Look at getLevelSetValues_SPH_singleThread for the actual implementation.
        levelSetValues.resize(getVertexCount(), -c);

        if (maxAllowedThreadSplits <= 0 || fluidParticles.size() <= maxAllowedThreadSplits) {
            getLevelSetValues_SPH_singleThread(&fluidParticles, 0, fluidParticles.size(), &neighborhood, h, progress);
            return;
        }
        
        std::vector<std::shared_ptr<std::thread>> threads;
        const unsigned long fluidParticlesPerThread = fluidParticles.size() / maxAllowedThreadSplits;

        for (int threadIndex = 0; threadIndex < maxAllowedThreadSplits; threadIndex++) {
            const uint64_t startIndex = fluidParticlesPerThread * threadIndex;
            const uint64_t endIndex = threadIndex == maxAllowedThreadSplits - 1 ? fluidParticles.size() : startIndex + fluidParticlesPerThread;

            std::shared_ptr<std::thread> thread = std::make_shared<std::thread>(&Voxelgrid::getLevelSetValues_SPH_singleThread, this, &fluidParticles, startIndex, endIndex, &neighborhood, h, progress);
            threads.emplace_back(thread);
        }

        for (auto& thread : threads) {
            thread->join();
        }
    }

private:
    void getLevelSetValues_SPH_singleThread(const std::vector<Eigen::Vector3d>* const fluidParticles, const uint64_t startIndex, const uint64_t endIndex, const learnSPH::Neighborhood* const neighborhood, const double h, ProgressData* const progress) {
        const double compactSupport = 2*h;

        for (int particleIndex = startIndex; particleIndex < endIndex; particleIndex++) {
            progress->particleIndex = particleIndex + 1; // +1 because it is zero based which looks weird in the output

            const Eigen::Vector3d& particle = (*fluidParticles)[particleIndex];

            // Check AABB
            std::vector<uint64_t> gridVerticesInCompactSupport;
            const Eigen::Vector3d lowerLeftBound = particle - Eigen::Vector3d(compactSupport, compactSupport, compactSupport);
            const Eigen::Vector3d upperRightBound = particle + Eigen::Vector3d(compactSupport, compactSupport, compactSupport);

            const std::array<uint64_t, 3> lowerLeftBoundCellCoordinates = { (uint64_t)((lowerLeftBound.x() - offset.x()) / cellEdgeLength), (uint64_t)((lowerLeftBound.y() - offset.y()) / cellEdgeLength), (uint64_t)((lowerLeftBound.z() - offset.z()) / cellEdgeLength) };
            const std::array<uint64_t, 3> upperRightBoundCellCoordinates = { (uint64_t)((upperRightBound.x() - offset.x()) / cellEdgeLength), (uint64_t)((upperRightBound.y() - offset.y()) / cellEdgeLength), (uint64_t)((upperRightBound.z() - offset.z()) / cellEdgeLength) };

            for (uint64_t i = lowerLeftBoundCellCoordinates.at(0); i < upperRightBoundCellCoordinates.at(0); i++) {
                for (uint64_t j = lowerLeftBoundCellCoordinates.at(1); j < upperRightBoundCellCoordinates.at(1); j++) {
                    for (uint64_t k = lowerLeftBoundCellCoordinates.at(2); k < upperRightBoundCellCoordinates.at(2); k++) {
                        const uint64_t vertexIndex = getVertexIndex(i, j, k);
                        const Eigen::Vector3d vertex = getVertexCoordinates(i, j, k);
                        
                        // Filter out grid vertices that are not in the compact support
                        const double distance = (vertex - particle).norm();
                        if (distance <= compactSupport) {
                            gridVerticesInCompactSupport.push_back(vertexIndex);
                        }
                    }
                }
            }

            // Add contribution of fluid particle to found grid vertices
            for (const uint64_t gridVertex : gridVerticesInCompactSupport) {
                uint64_t i, j, k;
                std::tie(i, j, k) = getijkFromV(gridVertex);

                double pi = 0;
                neighborhood->forEachNeighborDo(neighborhood->fluidPointSetId, particleIndex, neighborhood->fluidPointSetId, [&](const unsigned int neighborIndex) -> void {
                    pi += learnSPH::kernel::kernelCubicFunction((*fluidParticles)[particleIndex], (*fluidParticles)[neighborIndex]);
                });
                const double oneOverPi = 1.0 / pi;
                const double contribution = oneOverPi * learnSPH::kernel::kernelCubicFunction(getVertexCoordinates(i, j, k), (*fluidParticles)[particleIndex]);
                levelSetValues[gridVertex] += contribution;
            }
        }
    }

    inline GridCell getGridCell(uint64_t i, uint64_t j, uint64_t k) const {
        GridCell gridCell;
        for (int vertex = 0; vertex < CELL_VERTICES.size(); vertex++) {
            const Eigen::Vector3i vertexIndices(i + CELL_VERTICES[vertex][0], j + CELL_VERTICES[vertex][1], k + CELL_VERTICES[vertex][2]);
            const uint64_t V = getVertexIndex(vertexIndices.x(), vertexIndices.y(), vertexIndices.z());
            gridCell[vertex] = V;
        }
        return gridCell;
    }

    static std::tuple<uint64_t, uint64_t, uint64_t> getCellCount(const double boundingBoxWidth, const double boundingBoxHeight, const double boundingBoxDepth, const double cellEdgeLength) {
        const uint64_t cellsX = std::ceil(boundingBoxWidth / cellEdgeLength);
        const uint64_t cellsY = std::ceil(boundingBoxHeight / cellEdgeLength);
        const uint64_t cellsZ = std::ceil(boundingBoxDepth / cellEdgeLength);

        return std::tuple<uint64_t, uint64_t, uint64_t>(cellsX, cellsY, cellsZ);
    }
};
