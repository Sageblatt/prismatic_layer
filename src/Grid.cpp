#include "Grid.h"

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkWedge.h>

#include "Utils.h"

using std::cerr;
using std::endl;

using pl::Grid;
using pl::index_t;
using pl::get_normal_tr;

Grid::Grid(const std::string& filename, short int normal_sign):
        normal_sign(normal_sign) {
    vtkSmartPointer<vtkUnstructuredGrid> vtk_grid;
    std::string extension;
    if (filename.find_last_of('.') != std::string::npos)
        extension = filename.substr(filename.find_last_of('.'));

    // Drop the case of the extension
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   ::tolower);

    if (extension == ".vtu") {
        vtkNew<vtkXMLUnstructuredGridReader> reader;
        reader->SetFileName(filename.c_str());
        reader->Update();
        vtk_grid = reader->GetOutput();
    } else if (extension == ".vtk") {
        vtkNew<vtkUnstructuredGridReader> reader;
        reader->SetFileName(filename.c_str());
        reader->Update();
        vtk_grid = reader->GetOutput();
    } else {
        throw std::invalid_argument("Error: invalid extension of input file!");
    }

    auto vtk_faces = vtk_grid->GetCells();
    n_faces = vtk_grid->GetNumberOfCells();
    faces.resize(n_faces);

    for (index_t i = 0; i < n_faces; i++) {
        auto idarr = vtkSmartPointer<vtkIdList>::New();
        vtk_faces->GetCellAtId(i, idarr);
        for (auto j = 0; j < 3; j++) // TRIANGLES ARE HARDCODED!!!!
            faces[i](j) = idarr->GetId(j);
    }

    vtkSmartPointer<vtkPoints> vtk_pts = vtk_grid->GetPoints();
    n_pts = vtk_pts->GetNumberOfPoints();
    initial_pc.resize(n_pts);

    for (index_t i = 0; i < n_pts; i++) {
        auto p = vtk_pts->GetPoint(i);
        for (auto j = 0; j < 3; j++)
            initial_pc[i](j) = p[j];
    }

    connectivity.resize(initial_pc.size());

    for (index_t i = 0; i < n_faces; i++) {
        for (const auto& face_point : faces[i]) {
            for (const auto& neighbour : faces[i]) {
                if (neighbour == face_point)
                    continue;
                connectivity[face_point].push_back(neighbour);
            }
        }
    }
}

pair<vector<Vector3d>, vector<list<Matrix<double, 3, 1>>>>
Grid::getAllNormals(const vector<Vector3d>& point_cloud) const {
    auto sub_normals = vector<list<Matrix<double, 3, 1>>>(n_pts);

    for (index_t i = 0; i < n_faces; i++) {
        auto tri_pts = faces[i];
        auto tri_normal = get_normal_tr(point_cloud[tri_pts(0)],
                                        point_cloud[tri_pts(1)],
                                        point_cloud[tri_pts(2)],
                                        normal_sign);
        for (const auto point : faces[i])
            sub_normals[point].push_back(tri_normal);
    }

    auto normals = vector<Vector3d>(n_pts);
    for (index_t i = 0; i < n_pts; i++) {
        normals[i] = {0.0, 0.0, 0.0};
        for (const auto& sub_normal : sub_normals[i])
            normals[i] += sub_normal;
        normals[i] = normals[i] / normals[i].norm();
    }
    return {normals, sub_normals};
}

vector<Vector3d> Grid::getNormals(const vector<Vector3d>& point_cloud) const {
    return getAllNormals(point_cloud).first;
}

vector<Vector3d> Grid::laplacePointsSmooth(const vector<Vector3d>& points,
                                           const vector<Matrix<double, 3, 1>>& normals,
                                           double tau,
                                           unsigned iter,
                                           double H) {
    assert(tau > 0);
    assert(H > 0);

    auto smoothed = points;

    for (unsigned i = 0; i < iter; i++) {
        auto prev_iter_pts = smoothed;
        for (index_t j = 0; j < n_pts; j++) {
            auto neighbours = connectivity[j];
            auto current_point = smoothed[j];
            double total_weight = 0.0;
            Vector3d laplace = {0.0, 0.0, 0.0};
            for (auto const & neighbour : neighbours) {
                auto diff = current_point - prev_iter_pts[neighbour];
                double dst = diff.norm();
                total_weight += 1.0/dst;
                laplace += prev_iter_pts[neighbour] / dst;
            }
            smoothed[j] = prev_iter_pts[j] + tau *
                                             (1.0 / total_weight * laplace - prev_iter_pts[j]);
        }

        auto max_diff = H/iter;

        for (index_t j = 0; j < n_pts; j++) {
            auto diff = smoothed[j] - prev_iter_pts[j];
            auto product = diff.dot(normals[j]);
            if (product < 0) {
                smoothed[j] = prev_iter_pts[j];
            } else {
                if (product > max_diff)
                    product = max_diff;
                smoothed[j] = prev_iter_pts[j] + product * diff / diff.norm();
            }
        }
    }
    return smoothed;
}

void Grid::constructPL(unsigned layers_amount,
                       double multiplier,
                       double base,
                       unsigned iters_amount,
                       double tau) {
    assert(multiplier > 0);
    assert(base > 0);
    assert(tau > 0);
    assert(layers_amount > 2);

    processed_pc.clear();
    processed_pc.reserve(n_pts * layers_amount);

    auto normals = getNormals(initial_pc);
    auto new_layer = initial_pc;
    double Tm = 0.0;

    for (auto j = 0; j < layers_amount; j++) {
        for (index_t i = 0; i < n_pts; i++)
            new_layer[i] = new_layer[i] + Tm * normals[i];

        if (j != 0)
            new_layer = laplacePointsSmooth(new_layer, normals, tau, iters_amount, Tm);

        processed_pc.insert(processed_pc.end(), new_layer.begin(), new_layer.end());
        Tm = pow(1 + base, j) * multiplier;
        normals = getNormals(new_layer);
        // TODO: smoothing normals
    }
}

void Grid::exportPLToVTK(const string& filename) {
    auto vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto vtk_pts = vtkSmartPointer<vtkPoints>::New();
    vtk_pts->Allocate(processed_pc.size());
    for (auto const & pt : processed_pc)
        vtk_pts->InsertNextPoint(pt[0], pt[1], pt[2]);

    vtk_grid->SetPoints(vtk_pts);

    auto cells_size = n_faces; // Amount of cells in 1 layer
    auto offset = n_pts; // Amount of points in 1 layer
    unsigned n_layers = processed_pc.size() / n_pts;
    index_t n_cells = cells_size * (n_layers - 1);

    vtk_grid->AllocateExact(n_cells,n_cells * 6);

    auto prism = vtkSmartPointer<vtkWedge>::New();
    for (index_t i = 0; i < n_cells; i++) {
        auto layer = i / cells_size;
        prism->GetPointIds()->SetId(0,
                                    faces[i - layer*cells_size][0] + layer * offset);
        prism->GetPointIds()->SetId(1,
                                    faces[i - layer*cells_size][1] + layer * offset);
        prism->GetPointIds()->SetId(2,
                                    faces[i - layer*cells_size][2] + layer * offset);
        prism->GetPointIds()->SetId(3,
                                    faces[i - layer*cells_size][0] + (layer + 1) * offset);
        prism->GetPointIds()->SetId(4,
                                    faces[i - layer*cells_size][1] + (layer + 1) * offset);
        prism->GetPointIds()->SetId(5,
                                    faces[i - layer*cells_size][2] + (layer + 1) * offset);
        vtk_grid->InsertNextCell(prism->GetCellType(), prism->GetPointIds());
    }

    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(vtk_grid);
    writer->Write();
}

void Grid::exportNormalsToVTK(const string& filename) {
    auto normals = getNormals(initial_pc);
    auto vtk_normals = vtkSmartPointer<vtkDoubleArray>::New();
    vtk_normals->SetName("Normals");
    vtk_normals->SetNumberOfComponents(3);

    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    vtk_points->Allocate(initial_pc.size());

    for (index_t i = 0; i < n_pts; i++) {
        vtk_points->InsertNextPoint(initial_pc[i][0],
                                    initial_pc[i][1],
                                    initial_pc[i][2]);
        double arg[3] = {normals[i][0], normals[i][1], normals[i][2]};
        vtk_normals->InsertNextTuple(arg);
    }

    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->SetPoints(vtk_points);
    ug->GetPointData()->AddArray(vtk_normals);

    auto triangle = vtkSmartPointer<vtkTriangle>::New();
    for(index_t i = 0; i < n_faces; i++) {
        triangle->GetPointIds()->SetId(0, faces[i][0]);
        triangle->GetPointIds()->SetId(1, faces[i][1]);
        triangle->GetPointIds()->SetId(2, faces[i][2]);
        ug->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
    }

    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(ug);
    writer->Write();
}

