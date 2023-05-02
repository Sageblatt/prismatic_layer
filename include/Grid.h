#ifndef PRISMATIC_LAYER_GRID_H
#define PRISMATIC_LAYER_GRID_H

#include <vector>
#include <string>
#include <list>

#include <Eigen/Dense>

using std::vector;
using std::string;
using std::list;
using std::pair;

using Eigen::Vector3d;
using Eigen::Matrix;

namespace pl {
using index_t = unsigned long;

class Grid {
private:
    vector<Vector3d> initial_pc; // pc = point cloud
    index_t n_pts;

    vector<Matrix<index_t, 3, 1>> faces;
    index_t n_faces;

    vector<list<index_t>> connectivity;

    short int normal_sign;

    vector<Vector3d> processed_pc = vector<Vector3d>(0);

    pair<vector<Vector3d>, vector<list<Matrix<double, 3, 1>>>> getAllNormals(
            const vector<Vector3d>& point_cloud) const;
    vector<Vector3d> getNormals(const vector<Vector3d>& point_cloud) const;

    vector<Vector3d> laplacePointsSmooth(const vector<Vector3d>& points,
            const vector<Matrix<double, 3, 1>>& normals,
            double tau,
            unsigned iter,
            double H) const;

    vector<Vector3d> laplaceNormalsSmooth(const vector<Vector3d>& normals,
                                          double tau, unsigned iter) const;

public:
    explicit Grid(const string& filename, short int normal_sign);

    void constructPL(unsigned layers_amount,
                     double multiplier,
                     double base,
                     unsigned iters_amount,
                     double tau); // PL = prismatic layer

    void exportPLToVTK(const string& filename);

    void exportNormalsToVTK(const string& filename);
};

} // namespace pl

#endif //PRISMATIC_LAYER_GRID_H
