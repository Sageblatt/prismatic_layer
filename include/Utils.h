#ifndef PRISMATIC_LAYER_UTILS_H
#define PRISMATIC_LAYER_UTILS_H

#include <Eigen/Dense>

using Eigen::Vector3d;

namespace pl {

Vector3d get_normal_tr(const Vector3d& p1, const Vector3d& p2,
                              const Vector3d& p3, int sign);
} // namespace pl

#endif //PRISMATIC_LAYER_UTILS_H
