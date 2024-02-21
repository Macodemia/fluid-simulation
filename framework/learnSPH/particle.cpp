#include "particle.h"
#include <fstream>
#include "kernel.h"
#include "util/point_on_triangle.hpp"
#include "../extern/tiny_obj_loader/tiny_obj_loader.h"

std::tuple<std::vector<Eigen::Vector3d>, double, double> learnSPH::particle::sampleParticles(const double width, const double height, const double depth, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset, const Eigen::Vector3d& rotation) {
    std::vector<Eigen::Vector3d> box = {
            { 0, height, 0 },      // p_0
            { width, height, 0 },    // p_1
            { width, height, depth },   // p_2
            { 0, height, depth },    // p_3
            { 0, 0, 0 },        // p_4
            { width, 0, 0 },       // p_5
            { width, 0, depth },     // p_6
            { 0, 0, depth },      // p_7
    };

    return learnSPH::particle::sampleParticles(box, fluidRestDensity, numberOfParticles, offset, rotation);
}

// boxPositions is counter-clockwise. Starting at the top, then left. First point position is bottem left, looking from the front at the box.
// p_0: top area, front left
// p_1: top area, front right
// p_2: top area, back right
// p_3: top area, back left
// p_4: bottom area, front left
// p_5: bottom area, front right
// p_6: bottom area, back right
// p_7: bottom area, back left
std::tuple<std::vector<Eigen::Vector3d>, double, double> learnSPH::particle::sampleParticles(const std::vector<Eigen::Vector3d>& boxPositions,
                                                                                     const double fluidRestDensity,
                                                                                     const int numberOfParticles, const Eigen::Vector3d& offset, const Eigen::Vector3d& rotation) {
    if (numberOfParticles == 0) {
        return std::tuple<std::vector<Eigen::Vector3d>, double, double>(std::vector<Eigen::Vector3d>(), 0, 0);
    }

    const double width = (boxPositions[5] - boxPositions[4]).norm();
    const double height = (boxPositions[1] - boxPositions[5]).norm();
    const double depth = (boxPositions[6] - boxPositions[5]).norm();

    double volumeCuboid = width * height * depth;
    double massPerParticle = (volumeCuboid * fluidRestDensity) / numberOfParticles;
    double samplingDistance = std::pow(massPerParticle / fluidRestDensity, 1.0/3.0);
    
    std::vector<Eigen::Vector3d> particles;
    particles.reserve(numberOfParticles);
    
    for (double x = 0; x < width; x += samplingDistance) {
        for (double y = 0; y < height; y += samplingDistance) {
            for (double z = 0; z < depth; z += samplingDistance) {
                particles.emplace_back(x, y, z);
	    }
        }
    }

    if (numberOfParticles != particles.size()) {
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
       
        volumeCuboid = (maxX - minX + samplingDistance) * (maxY - minY + samplingDistance) * (maxZ - minZ + samplingDistance);
        massPerParticle = (volumeCuboid * fluidRestDensity) / particles.size();
        samplingDistance = std::pow(massPerParticle / fluidRestDensity, 1.0/3.0);
        
        std::cout << std::endl;
        std::cout << "Volume: " << volumeCuboid << std::endl;
        std::cout << "Number Of Particles: " << numberOfParticles << std::endl << "Actual Particle Count: " << particles.size() << std::endl;
    }

    for (Eigen::Vector3d& particle : particles) {
	if (rotation.x() != 0) {
	    particle = Eigen::AngleAxis<double>(rotation.x(), Eigen::Vector3d(1, 0, 0)) * particle;
	}
	if (rotation.y() != 0) {
	    //particle = particle * Eigen::AngleAxis<double>(rotation.x(), {1, 0, 0});
	}
	if (rotation.z() != 0) {
	    //particle = particle * Eigen::AngleAxis<double>(rotation.x(), {1, 0, 0});
	}
	particle = { particle.x() + offset.x(), particle.y() + offset.y(), particle.z() + offset.z() };
    }

    return std::tuple<std::vector<Eigen::Vector3d>, double, double>(particles, massPerParticle, samplingDistance / 2.0);
}

void learnSPH::particle::sampleParticles(learnSPH::PhysicalData& physicalData, const std::vector<Eigen::Vector3d>& boxPositions, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset, const Eigen::Vector3d& rotation) {
    std::tie(physicalData.fluidPositions, physicalData.fluidParticleMass, physicalData.particleRadius) = learnSPH::particle::sampleParticles(boxPositions, fluidRestDensity, numberOfParticles, offset, rotation);
    physicalData.velocities.resize(physicalData.fluidPositions.size(), Eigen::Vector3d(0, 0, 0));
    physicalData.densities.resize(physicalData.fluidPositions.size(), 0);
    physicalData.fluidSurfaceNormals.resize(physicalData.fluidPositions.size(), Eigen::Vector3d(0, 0, 0));
}

void learnSPH::particle::sampleParticles(learnSPH::PhysicalData& physicalData, const double width, const double height, const double depth, const double fluidRestDensity, const int numberOfParticles, const Eigen::Vector3d& offset, const Eigen::Vector3d& rotation) {
    std::tie(physicalData.fluidPositions, physicalData.fluidParticleMass, physicalData.particleRadius) = learnSPH::particle::sampleParticles(width, height, depth, fluidRestDensity, numberOfParticles, offset, rotation);
    physicalData.velocities.resize(physicalData.fluidPositions.size(), Eigen::Vector3d(0, 0, 0));
    physicalData.densities.resize(physicalData.fluidPositions.size(), 0);
    physicalData.fluidSurfaceNormals.resize(physicalData.fluidPositions.size(), Eigen::Vector3d(0, 0, 0));
}

std::vector<Eigen::Vector3d> learnSPH::particle::sampleTriangle(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C, const double samplingDistance, const bool hexagonalSampling) {
    // Handle obtuse triangles
    // by finding the longest edge and using its end points as A and B
    const double ABLength = (A - B).norm();
    const double ACLength = (A - C).norm();
    const double BCLength = (B - C).norm();
    if (ACLength > ABLength && ACLength > BCLength) {
        return sampleTriangle(A, C, B, samplingDistance, hexagonalSampling);
    } else if (BCLength > ABLength && BCLength > ACLength) {
        return sampleTriangle(B, C, A, samplingDistance, hexagonalSampling);
    }

    // We have correctly chosen A,B and C => Do sampling
    // i - Find Orthogonal vector basis on triangle plane
    const Eigen::Vector3d BMinusA = B - A;
    const Eigen::Vector3d CMinusA = C - A;
    const double halfSamplingDistance = samplingDistance / 2.0;

    const Eigen::Vector3d u = BMinusA.normalized();
    const Eigen::Vector3d triangleNormal = u.cross((CMinusA / CMinusA.norm())).normalized();
    const Eigen::Vector3d v = triangleNormal.cross(u).normalized();

    // ii - Find enlarging Vertices A', B', C'
    Eigen::Vector3d APrime = A;
    Eigen::Vector3d BPrime = B;
    Eigen::Vector3d CPrime = C;

    const Eigen::Vector3d triangleCenter = (A + B + C) / 3.0;
    APrime -= triangleCenter;
    BPrime -= triangleCenter;
    CPrime -= triangleCenter;

    Eigen::Matrix3d scalingMatrix;
    scalingMatrix <<
            1.0 + samplingDistance, 0, 0,
            0, 1.0 + samplingDistance, 0,
            0, 0, 1.0 + samplingDistance;

    APrime = scalingMatrix * APrime + triangleCenter;
    BPrime = scalingMatrix * BPrime + triangleCenter;
    CPrime = scalingMatrix * CPrime + triangleCenter;

    const Eigen::Vector3d CenterPrime = (APrime + BPrime + CPrime) / 3.0;

    // iii - Generate a collection of particles located at the nodes of the grid
    std::vector<Eigen::Vector3d> particles;

    if (hexagonalSampling) {
        int i = 0;
        const double xSamplingDistance = samplingDistance * std::sqrt(3.0) / 2.0;
        const double ySamplingDistance = samplingDistance / 2.0;
        const double ABSamplingSteps = std::ceil((APrime - BPrime).norm()) / xSamplingDistance;

        while (i < ABSamplingSteps) {
            int j = 0;
            const Eigen::Vector3d startPos = APrime + ySamplingDistance * v * (i % 2);

            while (true) {
                Eigen::Vector3d particle = startPos + u * xSamplingDistance * i + v * samplingDistance * j;

                if (!learnSPH::pointOnTriangle(APrime, BPrime, CPrime, triangleNormal, particle)) {
                    particle -= (particle - CenterPrime).normalized() * 1e-10;
                    if (!learnSPH::pointOnTriangle(APrime, BPrime, CPrime, triangleNormal, particle)) {
                        break;
                    }
                }

                particles.push_back(particle);
                j++;
            }
            i++;
        }

        return particles;
    } else {
        int i = 0;
        const double ABSamplingSteps = std::ceil((APrime - BPrime).norm()) / samplingDistance;

        while (i < ABSamplingSteps) {
            int j = 0;
            while (true) {
                Eigen::Vector3d particle = APrime + u * samplingDistance * i + v * samplingDistance * j;

                if (!learnSPH::pointOnTriangle(APrime, BPrime, CPrime, triangleNormal, particle)) {
                    particle -= (particle - CenterPrime).normalized() * 1e-10;
                    if (!learnSPH::pointOnTriangle(APrime, BPrime, CPrime, triangleNormal, particle)) {
                        break;
                    }
                }

                particles.push_back(particle);
                j++;
            }
            i++;
        }

        return particles;
    }
}

std::vector<Eigen::Vector3d> learnSPH::particle::sampleBoundaryBox(const double width, const double height, const double depth, const double samplingDistance, Eigen::Vector3d offset) {
    std::vector<Eigen::Vector3d> vertices = {
            { 0, height, 0 },      // p_0
            { width, height, 0 },    // p_1
            { width, height, depth },   // p_2
            { 0, height, depth },    // p_3
            { 0, 0, 0 },        // p_4
            { width, 0, 0 },       // p_5
            { width, 0, depth },     // p_6
            { 0, 0, depth },      // p_7
    };

    // Add offset
    for (Eigen::Vector3d& vertex : vertices) {
        vertex += offset;
    }

    // Clockwise
    std::vector<std::vector<int>> faces = {
            { 0, 1, 4 }, { 5, 4, 1 }, { 1, 5, 6 }, { 2, 6, 1 }, { 2, 6, 7 }, { 2, 3, 7 }, { 3, 0, 7 }, { 0, 4, 7 }, { 4, 7, 5 }, { 6, 7, 5 }
    };

    return learnSPH::particle::sampleBoundaryMesh(vertices, faces, samplingDistance);
}

std::vector<Eigen::Vector3d> learnSPH::particle::sampleBoundaryMesh(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::vector<int>>& faces, const double samplingDistance) {
    std::vector<Eigen::Vector3d> particles;

    for (std::vector<int> face : faces) {
        const std::vector<Eigen::Vector3d> triangle = learnSPH::particle::sampleTriangle(vertices[face[0]], vertices[face[1]], vertices[face[2]], samplingDistance);
        particles.insert(particles.end(), triangle.begin(), triangle.end());
    }

    return particles;
}

static void SplitString(std::string s, std::vector<std::string> &v){
    std::string temp = "";
    for(int i=0;i<s.length();++i){
        if(s[i]==' '){
            v.push_back(temp);
            temp = "";
        }
        else{
            temp.push_back(s[i]);
        }
    }
    v.push_back(temp);
}

std::vector<Eigen::Vector3d> learnSPH::particle::sampleObjFile(const std::string& file, const double samplingDistance, const Eigen::Vector3d& offset, const double scale) {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;

    std::ifstream stream;
    stream.open(file);
    //std::string objFileContent = "";
    std::string line;
    /*while(getline(stream, line)) {
        objFileContent += line + '\n';
    }*/
    std::vector<std::string> splittedLine;
#if 0
    tinyobj::ObjReaderConfig readerConfig;
    readerConfig.vertex_color = false;
    tinyobj::ObjReader reader;
    std::cout << "File: " << file << std::endl;
    reader.ParseFromString(objFileContent, "", readerConfig);
    std::cout << "Error: " << reader.Error() << std::endl;
    std::cout << "Warning: " << reader.Warning() << std::endl;

    vertices.resize(reader.GetAttrib().vertices.size());
    std::cout << "Vertices: " << reader.GetAttrib().vertices.size() << std::endl;
    for (unsigned int i = 0; i < reader.GetAttrib().vertices.size(); i++) {
        const tinyobj::real_t value = reader.GetAttrib().vertices[i];
        Eigen::Vector3d& vertex = vertices[i % 3];

        if (i % 3 == 0) {
            vertex.x() = value;    
        } else if (i % 3 == 1) {
            vertex.y() = value;
        } else if (i % 3 == 2) {
            vertex.z() = value;
        }
    }

    std::vector<tinyobj::shape_t> shapes = reader.GetShapes();
    std::cout << "Shapes: " << shapes.size() << std::endl;
    for (const tinyobj::shape_t& shape : shapes) {
        std::cout << "Shape-Faces: " << shape.mesh.indices.size() / 3.0 << std::endl;
        for (unsigned int i = 0; i < shape.mesh.indices.size(); i += 3) {
            const tinyobj::index_t& indexX = shape.mesh.indices[i];
            const tinyobj::index_t& indexY = shape.mesh.indices[i + 1];
            const tinyobj::index_t& indexZ = shape.mesh.indices[i + 2];
            faces.push_back({indexX.vertex_index, indexY.vertex_index, indexZ.vertex_index});
        }
    }
#endif

    while (getline(stream, line)) {
        if (line.length() == 0) continue;
        splittedLine.clear();

        if (line.at(0) == 'v') {
            SplitString(line, splittedLine);
            vertices.emplace_back(std::stod(splittedLine[1]) * scale + offset.x(), std::stod(splittedLine[2]) * scale + offset.y(), std::stod(splittedLine[3]) * scale + offset.z());
        } else if (line.at(0) == 'f') {
            SplitString(line, splittedLine);
            const int vIndex0 = std::stoi(splittedLine[1].substr(0, splittedLine[1].find('/'))) - 1;
            const int vIndex1 = std::stoi(splittedLine[2].substr(0, splittedLine[2].find('/'))) - 1;
            const int vIndex2 = std::stoi(splittedLine[3].substr(0, splittedLine[3].find('/'))) - 1;
            
            faces.push_back({vIndex0, vIndex1, vIndex2});
        }
    }
    stream.close();

    return learnSPH::particle::sampleBoundaryMesh(vertices, faces, samplingDistance);
}

