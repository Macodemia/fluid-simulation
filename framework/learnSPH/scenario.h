#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>

#include <Eigen/Dense>

enum SOLVER_TYPE { EXPLICIT, IMPLICIT };

// SOURCE: http://www.martinbroadhurst.com/how-to-trim-a-stdstring.html
static std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ") {
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
// SOURCE: http://www.martinbroadhurst.com/how-to-trim-a-stdstring.html
static std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ") {
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}
 
// SOURCE: http://www.martinbroadhurst.com/how-to-trim-a-stdstring.html
static std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ") {
    return ltrim(rtrim(str, chars), chars);
}

struct FluidBoxParameter {
    double width = 1;
    double height = 1;
    double depth = 1;
    int numberOfParticles = 1000;
    Eigen::Vector3d offset = Eigen::Vector3d(0.275, 2.0, 0.25);
    Eigen::Vector3d rotation = Eigen::Vector3d(0, 0, 0);
};

struct BoundaryBoxParameter {
    double width = 3.5;
    double height = 1.5;
    double depth = 1.5;
    double samplingDistance = 1/50.0;
    Eigen::Vector3d offset = Eigen::Vector3d(0, 0, 0);
};

struct BoundaryMeshParameter {
    std::string file;
    double samplingDistance = 1/50.0;
    Eigen::Vector3d offset = Eigen::Vector3d(0, 0, 0);
    double scale = 1.0;
};

constexpr bool printDebugInfo = false;

struct Scenario {
    SOLVER_TYPE solverType;
    std::string name = "Default Scenario";
    std::string filePath = "";
    double TOTAL_SECONDS_TO_SIMULATE = 1.0;
    double SAVE_FRAME_INTERVAL = 1.0/30.0; // ~33ms = ~30 fps
    double MAX_TIMESTEP_SIZE = 0.002; // 2ms
    std::vector<FluidBoxParameter> fluidBoxes;
    std::vector<BoundaryBoxParameter> boundaryBoxes;
    std::vector<BoundaryMeshParameter> boundaryMeshes;
    double FLUID_REST_DENSITY = 1000; // kg/m^3
    double ETA = 1.2; // tuning parameter for the smoothing length
    double EPSILON = 0.5; // Smoothed particle velocities
    double STIFFNESS_B = 1000;
    double VISCOSITY_V_f = 0.05; // viscosity value of the fluid
    double VISCOSITY_V_b = 0.0; // viscosity value of the boundary (friction)
    double GRAVITY = 9.81;
    double MAX_FLUID_PARTICLE_VELOCITY = 3.0; // m/s
    int CONSTRAINT_ITERATIONS = 4; // ONLY for implicit solver
    double VELOCITY_TIME_STEP_FACTOR = 0.5;
    double COHESION = 0.0;
    double ADHESION = 0.0;

    // TODO [low-prio]: This function is kind of a (big)mess
    void readFromFile(std::string scenarioFile) {
        std::ifstream stream;
        stream.open(scenarioFile);
        if (!stream.is_open()) {
            std::cout << "ERROR: readFromFile (scenario) failed! Passed file (" << scenarioFile << ") could not be opened!" << std::endl;
            return;
        } 
        
        if (printDebugInfo) {
            std::cout << "Parsing Scenario File: " << scenarioFile << std::endl;
        }

        const int lastSlash = scenarioFile.find_last_of('/');
        //this->name = std::filesystem::path(scenarioFile).stem();
        const std::string withoutDirectories = scenarioFile.substr(lastSlash + 1);
        const int lastDot = withoutDirectories.find_last_of('.');
        this->name = withoutDirectories.substr(0, lastDot);
        this->filePath = scenarioFile;

        // Clear default vectors
        fluidBoxes.clear();
        boundaryBoxes.clear();

        std::string line;
        bool readingFluidBox = false;
        bool readingBoundaryBox = false;
        bool readingBoundaryMesh = false;
        while(getline(stream, line)) {
            line = trim(line);
            if (line.empty()) {
                readingFluidBox = readingBoundaryBox = readingBoundaryMesh = false;
                if (printDebugInfo) {
                    std::cout << "INFO: Skipping empty line: " << line << std::endl;
                    std::cout << "INFO: Disabling FluidBox, BoundaryBox and BoundaryMesh reading mode!" << std::endl;
                }
                continue;
            }

            if (line.find("explicit") == 0 || line.find("EXPLICIT") == 0) {
                 if (solverType == SOLVER_TYPE::EXPLICIT) {
                    line = line.erase(0, std::string("explicit").length());
                 } else {
                    continue;
                 }
            }
            else if (line.find("implicit") == 0 || line.find("IMPLICIT") == 0) {
                 if (solverType == SOLVER_TYPE::IMPLICIT) {
                    line = line.erase(0, std::string("implicit").length());
                 } else {
                    continue;
                 }
            }
            line = trim(line);

            int indexOfColumn = line.find(':');
            if (indexOfColumn == std::string::npos) {
                std::cout << "WARNING: readFromFile (scenario) failed! Line (" << line << ") does not contain a ':'. File: " << scenarioFile << std::endl;
                continue;
            }

            std::string paramName = line.substr(0, indexOfColumn);
            trim(paramName);
            std::string paramValue = line.substr(indexOfColumn + 1);
            trim(paramValue);

            if (readingFluidBox) {
                FluidBoxParameter& fluidBox = fluidBoxes[fluidBoxes.size() - 1];
                if (paramName == "width") {
                    fluidBox.width = std::stod(paramValue);
                } else if (paramName == "height") {
                    fluidBox.height = std::stod(paramValue);
                } else if (paramName == "depth") {
                    fluidBox.depth = std::stod(paramValue);
                } else if (paramName == "numberOfParticles") {
                    fluidBox.numberOfParticles = std::stoi(paramValue);
                } else if (paramName == "offset") {
                    const int firstSpace = paramValue.find(' ');
                    const int secondSpace = paramValue.find(' ', firstSpace + 1);
                    const double offsetX = std::stod(paramValue.substr(0, firstSpace));
                    const double offsetY = std::stod(paramValue.substr(firstSpace + 1, secondSpace));
                    const double offsetZ = std::stod(paramValue.substr(secondSpace));
                    fluidBox.offset = Eigen::Vector3d(offsetX, offsetY, offsetZ);
                } else if (paramName == "rotation") {
                    const int firstSpace = paramValue.find(' ');
                    const int secondSpace = paramValue.find(' ', firstSpace + 1);
                    const double rotationX = std::stod(paramValue.substr(0, firstSpace));
                    const double rotationY = std::stod(paramValue.substr(firstSpace + 1, secondSpace));
                    const double rotationZ = std::stod(paramValue.substr(secondSpace));
                    fluidBox.rotation = Eigen::Vector3d(rotationX, rotationY, rotationZ);
		} else {
                    std::cout << "WARNING: readFromFile (scenario) failed! ParamName (" << paramName << ") is not a valid parameter for fluid box. File: " << scenarioFile << std::endl;
                }      
            } else if (readingBoundaryBox) {
                BoundaryBoxParameter& boundaryBox = boundaryBoxes[boundaryBoxes.size() - 1];

                if (paramName == "width") {
                    boundaryBox.width = std::stod(paramValue);
                } else if (paramName == "height") {
                    boundaryBox.height = std::stod(paramValue);
                } else if (paramName == "depth") {
                    boundaryBox.depth = std::stod(paramValue);
                } else if (paramName == "samplingDistance") {
                    boundaryBox.samplingDistance = std::stod(paramValue);
                } else if (paramName == "offset") {
                    const int firstSpace = paramValue.find(' ');
                    const int secondSpace = paramValue.find(' ', firstSpace + 1);
                    const double offsetX = std::stod(paramValue.substr(0, firstSpace));
                    const double offsetY = std::stod(paramValue.substr(firstSpace + 1, secondSpace));
                    const double offsetZ = std::stod(paramValue.substr(secondSpace));
                    boundaryBox.offset = Eigen::Vector3d(offsetX, offsetY, offsetZ);
                } else {
                    std::cout << "WARNING: readFromFile (scenario) failed! ParamName (" << paramName << ") is not a valid parameter for boundary box. File: " << scenarioFile << std::endl;
                }
            } else if (readingBoundaryMesh) {
                BoundaryMeshParameter& boundaryMesh = boundaryMeshes[boundaryMeshes.size() - 1];

                if (paramName == "file") {
                    boundaryMesh.file = paramValue;
                } else if (paramName == "samplingDistance") {
                    boundaryMesh.samplingDistance = std::stod(paramValue);
                } else if (paramName == "scale") {
                    boundaryMesh.scale = std::stod(paramValue);
                } else if (paramName == "offset") {
                    const int firstSpace = paramValue.find(' ');
                    const int secondSpace = paramValue.find(' ', firstSpace + 1);
                    const double offsetX = std::stod(paramValue.substr(0, firstSpace));
                    const double offsetY = std::stod(paramValue.substr(firstSpace + 1, secondSpace));
                    const double offsetZ = std::stod(paramValue.substr(secondSpace));
                    boundaryMesh.offset = Eigen::Vector3d(offsetX, offsetY, offsetZ);
                } else {
                    std::cout << "WARNING: readFromFile (scenario) failed! ParamName (" << paramName << ") is not a valid parameter for boundary mesh. File: " << scenarioFile << std::endl;
                }
            } else {
                if (paramName == "FluidBox") {
                    if (printDebugInfo) {
                        std::cout << "INFO: Enabling FluidBox reading mode!" << std::endl;
                    }
                    readingFluidBox = true;
                    fluidBoxes.push_back(FluidBoxParameter());
                    continue;
                } else if (paramName == "BoundaryBox") {
                    if (printDebugInfo) {
                        std::cout << "INFO: Enabling BoundaryBox reading mode!" << std::endl;
                    }
                    readingBoundaryBox = true;
                    boundaryBoxes.push_back(BoundaryBoxParameter());
                    continue;
                } else if (paramName == "BoundaryMesh") {
                    if (printDebugInfo) {
                        std::cout << "INFO: Enabling BoundaryMesh reading mode!" << std::endl;
                    }
                    readingBoundaryMesh = true;
                    boundaryMeshes.push_back(BoundaryMeshParameter());
                } else if (paramName == "TOTAL_SECONDS_TO_SIMULATE") {
                   TOTAL_SECONDS_TO_SIMULATE = std::stod(paramValue);  
                } else if (paramName == "SAVE_FRAME_INTERVAL") {
                    SAVE_FRAME_INTERVAL = std::stod(paramValue);
                } else if (paramName == "MAX_TIMESTEP_SIZE") {
                    MAX_TIMESTEP_SIZE = std::stod(paramValue);
                } else if (paramName == "FLUID_REST_DENSITY") {
                    FLUID_REST_DENSITY = std::stod(paramValue);
                } else if (paramName == "ETA") {
                    ETA = std::stod(paramValue);
                } else if (paramName == "EPSILON") {
                    EPSILON = std::stod(paramValue);
                } else if (paramName == "STIFFNESS_B") {
                    STIFFNESS_B = std::stod(paramValue);
                } else if (paramName == "VISCOSITY_V_f") {
                    VISCOSITY_V_f = std::stod(paramValue);
                } else if (paramName == "VISCOSITY_V_b") {
                    VISCOSITY_V_b = std::stod(paramValue);
                } else if (paramName == "GRAVITY") {
                    GRAVITY = std::stod(paramValue);
                } else if (paramName == "MAX_FLUID_PARTICLE_VELOCITY") {
                    MAX_FLUID_PARTICLE_VELOCITY = std::stod(paramValue);
                } else if (paramName == "CONSTRAINT_ITERATIONS") {
                    CONSTRAINT_ITERATIONS = std::stoi(paramValue);
                } else if (paramName == "VELOCITY_TIME_STEP_FACTOR") {
                    VELOCITY_TIME_STEP_FACTOR = std::stod(paramValue);
                } else if (paramName == "COHESION") {
                    COHESION = std::stod(paramValue);
                } else if (paramName == "ADHESION") {
                    ADHESION = std::stod(paramValue);
                } else {
                    std::cout << "WARNING: readFromFile (scenario) failed! ParamName (" << paramName << ") is not a valid parameter. File: " << scenarioFile << std::endl;
                }
            }
            
            if (printDebugInfo) {
                std::cout << "ParamName: " << paramName << std::endl;
                std::cout << "ParamValue: " << paramValue << std::endl;
            }
        }  

        stream.close();
    }

    std::string getSolverTypeString() const {
        switch (solverType) {
            case SOLVER_TYPE::EXPLICIT:
                return "explicit";
            case SOLVER_TYPE::IMPLICIT:
                return "implicit";
            default:
                return "UNRECOGNIZED_SOLVER";
        }
    }

    void print() const {
        std::cout << "--------------------------" << std::endl;
        std::cout << "Scenario: " << name << " (" << filePath << ")" << std::endl;

        std::cout << "\tTime Parameters: " << std::endl;
        std::cout << "\t\tTOTAL_SECONDS_TO_SIMULATE: " << TOTAL_SECONDS_TO_SIMULATE << std::endl;
        std::cout << "\t\tSAVE_FRAME_INTERVAL: " << SAVE_FRAME_INTERVAL << std::endl;
        std::cout << "\t\tMAX_TIMESTEP_SIZE: " << MAX_TIMESTEP_SIZE << std::endl;
        std::cout << "\t\tCONSTRAINT_ITERATIONS: " << CONSTRAINT_ITERATIONS << std::endl;
        std::cout << "\t\tVELOCITY_TIME_STEP_FACTOR: " << VELOCITY_TIME_STEP_FACTOR << std::endl;
        std::cout << std::endl;

        std::cout << "\tFluid Parameters: " << std::endl;
        std::cout << "\t\tFLUID_REST_DENSITY: " << FLUID_REST_DENSITY << std::endl;
        std::cout << "\t\tETA: " << ETA << std::endl;
        std::cout << "\t\tEPSILON: " << EPSILON << std::endl;
        std::cout << "\t\tSTIFFNESS_B: " << STIFFNESS_B << std::endl;
        std::cout << "\t\tVISCOSITY_V_f: " << VISCOSITY_V_f << std::endl;
        std::cout << "\t\tVISCOSITY_V_b: " << VISCOSITY_V_b << std::endl;
        std::cout << "\t\tGRAVITY: " << GRAVITY << std::endl;
        std::cout << "\t\tMAX_FLUID_PARTICLE_VELOCITY: " << MAX_FLUID_PARTICLE_VELOCITY << std::endl;
        std::cout << "\t\tCOHESION: " << COHESION << std::endl;
        std::cout << "\t\tADHESION: " << ADHESION << std::endl;

        std::cout << "\tFluid Boxes" << std::endl;
        for (int i = 0; i < fluidBoxes.size(); i++) {
            const FluidBoxParameter& fluidBox = fluidBoxes[i];
            std::cout << "\t\t" << i << ":" << std::endl;
            std::cout << "\t\twidth: " << fluidBox.width << std::endl;
            std::cout << "\t\theight: " << fluidBox.height << std::endl;
            std::cout << "\t\tdepth: " << fluidBox.depth << std::endl;
            std::cout << "\t\tnumberOfParticles: " << fluidBox.numberOfParticles << std::endl;
            std::cout << "\t\toffset: (" << fluidBox.offset.x() << "," << fluidBox.offset.y() << "," << fluidBox.offset.z() << ")" << std::endl;
            std::cout << "\t\trotation: (" << fluidBox.rotation.x() << "," << fluidBox.rotation.y() << "," << fluidBox.rotation.z() << ")" << std::endl;
        }

        std::cout << "\tBoundary Boxes" << std::endl;
        for (int i = 0; i < boundaryBoxes.size(); i++) {
            const BoundaryBoxParameter& boundaryBox = boundaryBoxes[i];
            std::cout << "\t\t" << i << ":" << std::endl;
            std::cout << "\t\twidth: " << boundaryBox.width << std::endl;
            std::cout << "\t\theight: " << boundaryBox.height << std::endl;
            std::cout << "\t\tdepth: " << boundaryBox.depth << std::endl;
            std::cout << "\t\tsamplingDistance: " << boundaryBox.samplingDistance << std::endl;
            std::cout << "\t\toffset: (" << boundaryBox.offset.x() << "," << boundaryBox.offset.y() << "," << boundaryBox.offset.z() << ")" << std::endl;
        }

        std::cout << "\tBoundary Meshes" << std::endl;
        for (int i = 0; i < boundaryMeshes.size(); i++) {
            const BoundaryMeshParameter& boundaryMesh = boundaryMeshes[i];
            std::cout << "\t\t" << i << ":" << std::endl;
            std::cout << "\t\tfile: " << boundaryMesh.file << std::endl;
            std::cout << "\t\tsamplingDistance: " << boundaryMesh.samplingDistance << std::endl;
            std::cout << "\t\toffset: (" << boundaryMesh.offset.x() << "," << boundaryMesh.offset.y() << "," << boundaryMesh.offset.z() << ")" << std::endl;
            std::cout << "\t\tscale: " << boundaryMesh.scale << std::endl;
        }

        std::cout << "--------------------------" << std::endl;
    }
};
