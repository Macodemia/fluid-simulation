#pragma once
#include <vector>

#include <Eigen/Dense>

namespace learnSPH
{
	static bool pointOnTriangle(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C, const Eigen::Vector3d& normal, const Eigen::Vector3d& point)
	{
		const Eigen::Vector3d pab = (B - A).cross(point - A);
		if (pab.dot(normal) < 0.0) { return false; }

		const Eigen::Vector3d pbc = (C - B).cross(point - B);
		if (pbc.dot(normal) < 0.0) { return false; }

		const Eigen::Vector3d pca = (A - C).cross(point - C);
		if (pca.dot(normal) < 0.0) { return false; }

		return true;
	}
}