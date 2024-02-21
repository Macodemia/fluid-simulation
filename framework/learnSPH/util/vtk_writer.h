#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <vtkio/VTKFile.h>

namespace learnSPH
{
	void save_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles, const std::vector<double>& particle_scalar_data, const std::vector<Eigen::Vector3d>& particle_vector_data);
	void saveTriMeshToVTK(std::string path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles);

	void _compute_node_normals(std::vector<Eigen::Vector3d>& out_normals, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles);
}
