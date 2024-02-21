#include <iostream>
#include <cstdio>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <functional>
#include <thread>
#include "voxelgrid.h"
#include "../learnSPH/util/vtk_writer.h"
#include "../learnSPH/kernel.h"
#include "../learnSPH/particle.h"
#include "../learnSPH/neighborhood.h"

/*
 * Example:
---------------------------------------------------------------- Progress ----------------------------------------------------------------
| Thread 0                        | Thread 1                        | Thread 2                        | Thread 3                        |
| Processing frameFile: 10/37     | Processing frameFile: 10/37     | Processing frameFile: 10/37     | Processing frameFile: 10/37     |
| Calc contribution: 5000/2100123 | Calc contribution: 5000/2100123 | Calc contribution: 5000/2100123 | Calc contribution: 5000/2100123 |
| Thread 4                        | Thread 5                        | Thread 6                        | Thread 7                        |
| Processing frameFile: 10/37     | Processing frameFile: 10/37     | Processing frameFile: 10/37     | Processing frameFile: 10/37     |
| Calc contribution: 5000/2100123 | Calc contribution: 5000/2100123 | Calc contribution: 5000/2100123 | Calc contribution: 5000/2100123 |
------------------------------------------------------------------------------------------------------------------------------------------
 */
static bool firstTimePrinted = true;
void printProgress(const std::vector<ProgressData>& progress) {
    constexpr int tableColumns = 4;
    constexpr int tableCellWidth = 30;
    const int tableRows = std::ceil((double)progress.size() / (double)tableColumns);

    const std::string titleDashes((tableColumns * tableCellWidth + tableColumns) / 2, '-');
    const std::string title = " Progress ";
    const std::string endDashes(titleDashes.length() * 2 + title.length(), '-');

    if (firstTimePrinted) {
        firstTimePrinted = false;
        // Clear entire screen and move cursor to top left
        printf("\e[2J");
        //printf("\e[s");
    } else {
        //printf("\e[u");
    }

    // Move Cursor to Top Left
    printf("\e[H");

    std::cout << titleDashes << title << titleDashes << std::endl;

    for (int row = 0; row < tableRows; row++) {
        // threadNumber row
        printf(" ");
        for (int i = row * tableColumns; i < std::min(row * tableColumns + tableColumns, (int)progress.size()); i++) {
            const ProgressData& progressData = progress[i];

            const std::string threadNumberString = "Thread " + std::to_string(i);
            printf("| %-*s|", tableCellWidth, threadNumberString.c_str());
        }

        std::cout << std::endl;

        // frameFile progress row
        printf(" ");
        for (int i = row * tableColumns; i < std::min(row * tableColumns + tableColumns, (int)progress.size()); i++) {
            const ProgressData& progressData = progress[i];

            const std::string frameFileString = "frameFile: " + progressData.frameFile + "/" + std::to_string(progressData.endIndex);
            printf("| %-*s|", tableCellWidth, frameFileString.c_str());
        }

        std::cout << std::endl;

        // particleIndex progress row
        printf(" ");
        for (int i = row * tableColumns; i < std::min(row * tableColumns + tableColumns, (int)progress.size()); i++) {
            const ProgressData& progressData = progress[i];

            const std::string particleIndexString = "particleIndex " + std::to_string(progressData.particleIndex) + "/" + std::to_string(progressData.particleCount);
            printf("| %-*s|", tableCellWidth, particleIndexString.c_str());
        }

        std::cout << std::endl << endDashes << std::endl;
    }
}

double capitalPhi_sphere(const Eigen::Vector3d& point, const double radius) {
    return point.norm() - radius;
}

double capitalPhi_torus(const Eigen::Vector3d& point, const double minorRadius, const double majorRadius) {
    return (minorRadius*minorRadius) - (sqrt(point.x() * point.x() + point.y() * point.y()) - majorRadius) * (sqrt(point.x() * point.x() + point.y() * point.y()) - majorRadius) - (point.z() * point.z());
}

void getLevelSetValues_Shape(const Voxelgrid& voxelgrid, std::vector<double>& levelSetValues, const std::function<double(const Eigen::Vector3d&)> capitalPhi) {
    for (int V = 0; V < voxelgrid.getVertexCount(); V++) {
        uint64_t i, j, k;
        std::tie(i, j, k) = voxelgrid.getijkFromV(V);
        levelSetValues[V] = capitalPhi(voxelgrid.getVertexCoordinates(i, j, k));
    }
}

void parseFluidParticleFile(const std::string& file, std::vector<Eigen::Vector3d>& fluidParticles, double& h, std::string& scenarioName) {
    std::ifstream frameFile;
    frameFile.open(file);
    if (!frameFile.is_open()) {
        std::cout << "Error: Could not load frame file. Filename: " << file;
        std::cout << "Aborting!" << std::endl;
        return;
    }

    std::string line;
    // Read scenarioName
    getline(frameFile, line);
    scenarioName = line;

    // Read h
    getline(frameFile, line);
    h = std::stod(line);

    // Read particle radius
    getline(frameFile, line);
    const double particleRadius = std::stod(line); // TODO: Probably not needed -> delete

    // Read particles
    while (getline(frameFile, line)) {
        const int firstSpace = line.find(' ');
        const int secondSpace = line.find(' ', firstSpace + 1);
        std::string xString = line.substr(0, firstSpace);
        std::string yString = line.substr(firstSpace + 1, secondSpace - firstSpace - 1);
        std::string zString = line.substr(secondSpace + 1, line.length() - secondSpace - 1);
        fluidParticles.emplace_back(std::stod(xString), std::stod(yString), std::stod(zString));
    }
    frameFile.close();
}

void computeMinimalBoundingBox(const std::vector<Eigen::Vector3d>& particles, const double h, Eigen::Vector3d& offset, double& boundingBoxWidth, double& boundingBoxHeight, double& boundingBoxDepth) {
    // Calculate smallest bounding box from particles
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

    const double safetyMargin = (2 * h);
    boundingBoxWidth = maxX - minX + safetyMargin * 2;
    boundingBoxHeight = maxY - minY + safetyMargin * 2;
    boundingBoxDepth = maxZ - minZ + safetyMargin * 2;
    offset = Eigen::Vector3d(minX - safetyMargin, minY - safetyMargin, minZ - safetyMargin);
}

void processFrameFile(const std::vector<std::string>& frameFiles, const int startIndex, const int endIndex, const double cellEdgeLength, const double C, ProgressData* const progress, int maxAllowedThreadSplits = 0) {
    for (int fileIndex = startIndex; fileIndex < endIndex; fileIndex++) {
        const std::string& file = frameFiles[fileIndex];
        const std::string stringToFind = "fluid_for_reconstruction_";
        progress->frameFile = file.substr(file.find(stringToFind) + stringToFind.length());
        progress->endIndex = endIndex;

        std::string scenarioName;
        std::vector<Eigen::Vector3d> fluidParticles;
        double h;
        parseFluidParticleFile(file, fluidParticles, h, scenarioName);
        progress->particleCount = fluidParticles.size();

        std::vector<Eigen::Vector3d> boundaryParticles;
        learnSPH::kernel::computeLookupTables(h);
        learnSPH::Neighborhood neighborhood(h, fluidParticles, boundaryParticles);

        Eigen::Vector3d boundingBoxOffset;
        double boundingBoxWidth, boundingBoxHeight, boundingBoxDepth;
        computeMinimalBoundingBox(fluidParticles, h, boundingBoxOffset, boundingBoxWidth, boundingBoxHeight, boundingBoxDepth);

        const std::string gridVerticesFileName = "grid_vertices_" + std::to_string(fileIndex) + ".vtk";
        Voxelgrid voxelgrid;
        voxelgrid.init(boundingBoxOffset, boundingBoxWidth, boundingBoxHeight, boundingBoxDepth, cellEdgeLength, gridVerticesFileName);
        voxelgrid.getLevelSetValues_SPH(fluidParticles, neighborhood, h, C, progress, maxAllowedThreadSplits);

        std::vector<Eigen::Vector3d> intersectionPoints;
        std::vector<std::array<int, 3>> triangles;
        std::tie(intersectionPoints, triangles) = voxelgrid.marchingCubes();
        
        if (!intersectionPoints.empty()) {
            const std::string basePath = "../res/" + scenarioName + "/reconstruction_" + progress->frameFile;//+ std::to_string(fileIndex);
            std::string vtkFilePath = basePath + ".vtk";
            learnSPH::saveTriMeshToVTK(vtkFilePath, intersectionPoints, triangles);
#if 0
            // save as obj
            std::ofstream objFileStream;
            objFileStream.open(basePath + ".obj");
            for (const auto& vertex : intersectionPoints) {
                objFileStream << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
            }

            for (const auto& face : triangles) {
                if (face[0] >= intersectionPoints.size() || face[1] >= intersectionPoints.size() || face[2] >= intersectionPoints.size()) {
                    std::cout << "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                    std::cout << "Intersection Point Size: " << intersectionPoints.size() << std::endl;
                    std::cout << "Face: " << face[0] << " " << face[1] << " " << face[2] << std::endl;
                }
                objFileStream << "f " << face[0]+1 << " " << face[1]+1 << " " << face[2]+1 << std::endl;
            }
            objFileStream.close();
#endif
        }
    }

    progress->threadFinished = true;
}

// Commandline Arguments:
// A list of frame files
// Example execution:
// ./learnSPH_reconstruction build/release/app/reconstruction_*
int main(int argc, char** argv) {
    const double cellEdgeLength = 0.01;
    const double C = 0.6;
   
    /* 
    //Reconstructing a Sphere/Torus 
    Voxelgrid voxelgrid;
    voxelgrid.init({-1,-1,-1},2,2,2,cellEdgeLength,"sphere");
    
    for (const auto& point : voxelgrid.gridVertices) {
        voxelgrid.levelSetValues.push_back(capitalPhi_sphere(point, 0.5));
        //voxelgrid.levelSetValues.push_back(capitalPhi_torus(point, 0.2, 0.5));
    }

    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int, 3>> indices;
    std::tie(vertices, indices) = voxelgrid.marchingCubes();
    learnSPH::saveTriMeshToVTK("sphere.vtk", vertices, indices);

    return 0;
    */

    if (argc == 1) {
        std::cout << "Please specify frame files as parameters!" << std::endl;
        return 0;
    }

    std::vector<std::string> frameFiles;
    // Skip first element, as it is the executable path
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        frameFiles.push_back(arg);
        std::cout << "Adding frameFile to be reconstructed: " << arg << std::endl;
    }

    const unsigned int numberOfThreads = 3;
    const unsigned int filesPerThread = frameFiles.size() / numberOfThreads;
    const unsigned int numberOfThreadSplits = 0;//(std::thread::hardware_concurrency() - numberOfThreads) / numberOfThreads;

    std::cout << "Number of Threads: " << numberOfThreads << std::endl;
    std::cout << "Number of ThreadSplits: " << numberOfThreadSplits << std::endl;
    std::vector<std::shared_ptr<std::thread>> threads;

    std::vector<ProgressData> progress(numberOfThreads);

    for (int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++) {
        const int startIndex = filesPerThread * threadIndex;
        const int endIndex = threadIndex == numberOfThreads - 1 ? frameFiles.size() : startIndex + filesPerThread;

        std::shared_ptr<std::thread> thread = std::make_shared<std::thread>(processFrameFile, frameFiles, startIndex, endIndex, cellEdgeLength, C, &progress[threadIndex], numberOfThreadSplits);
        threads.emplace_back(thread);
    }

    bool allThreadsFinished = false;
    while (!allThreadsFinished) {
        allThreadsFinished = true;
        for (const ProgressData& progressData : progress) {
            if (!progressData.threadFinished) {
                allThreadsFinished = false;
                break;
            }
        }

        printProgress(progress);
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    for (auto& thread : threads) {
        thread->join();
    }
    
    return 0;
}
