// STL Headers
#include <iostream>
#include <algorithm>
#include <string>
#include <list>
#include <vector>

// VTK Headers
#include <vtkAppendFilter.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkWedge.h>
#include <vtkDataWriter.h>
//#include <vtkFeatureEdges.h>

// Eigen3 Headers
#include <Eigen/Dense>


using std::cout;
using std::endl;
using std::list;
using std::vector;
using std::pair;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Vector3d;


vtkSmartPointer<vtkUnstructuredGrid> ReadUnstructuredGrid(
        std::string const& fileName) {
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  std::string extension = "";
  if (fileName.find_last_of(".") != std::string::npos)
  {
    extension = fileName.substr(fileName.find_last_of("."));
  }

  // Drop the case of the extension
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 ::tolower);

  if (extension == ".vtu")
  {
    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(fileName.c_str());
    reader->Update();
    unstructuredGrid = reader->GetOutput();
  }
  else if (extension == ".vtk")
  {
    vtkNew<vtkUnstructuredGridReader> reader;
    reader->SetFileName(fileName.c_str());
    reader->Update();
    unstructuredGrid = reader->GetOutput();
  }
  else
  {
    vtkNew<vtkSphereSource> source;
    source->Update();
    vtkNew<vtkAppendFilter> appendFilter;
    appendFilter->AddInputData(source->GetOutput());
    appendFilter->Update();
    unstructuredGrid = appendFilter->GetOutput();
  }

  return unstructuredGrid;
}

template <typename T>
void print(T& input) {
    cout << input << endl;
}

template <typename T>
void print(T&& input) {
    cout << input << endl;
}

Vector3d get_normal_tr(const Vector3d& p1,
                       const Vector3d& p2,
                       const Vector3d& p3,
                       int sign) {
    auto v1 = p3-p2;
    auto v2 = p1-p2;
    return v1.cross(v2) * sign;
}

pair<vector<Vector3d>, vector<list<Matrix<double, 3, 1>>>> get_all_normals(
        const vector<Vector3d>& point_cloud,
        const vector<Matrix<unsigned long, 3, 1>>& faces,
        int normal_sign) {
    auto n_pts = point_cloud.size();
    auto n_faces = faces.size();

    auto sub_normals = vector<list<Matrix<double, 3, 1>>>(n_pts);

    for(auto i = 0; i < n_faces; i++) {
        auto tri_pts = faces[i];
        auto tri_normal = get_normal_tr(point_cloud[tri_pts(0)],
                                        point_cloud[tri_pts(1)],
                                        point_cloud[tri_pts(2)],
                                        normal_sign);
        for(const auto point : faces[i])
            sub_normals[point].push_back(tri_normal);
    }

    auto normals = vector<Vector3d>(n_pts);
    for(auto i = 0; i < n_pts; i++) {
        normals[i] = {0.0, 0.0, 0.0};
        for(const auto& sub_normal : sub_normals[i])
            normals[i] += sub_normal;
        normals[i] = normals[i] / normals[i].norm();
    }
    return {normals, sub_normals};
}

vector<Vector3d> get_normals(const vector<Vector3d>& point_cloud,
                             const vector<Matrix<unsigned long, 3, 1>>& faces,
                             int normal_sign) {
    return get_all_normals(point_cloud, faces, normal_sign).first;
}

vector<list<unsigned long>> get_connectivity_array(
        unsigned long number_of_points,
        const vector<Matrix<unsigned long, 3, 1>>& faces) {
    vector<list<unsigned long>> result(number_of_points);
    auto n = faces.size();

    for(auto i = 0; i < n; i++) {
        for(const auto& face_point : faces[i]) {
            for(const auto& neighbour : faces[i]) {
                if(neighbour == face_point)
                    continue;
                result[face_point].push_back(neighbour);
            }
        }
    }
    return result;
}

vector<Vector3d> laplace_points_smooth(const vector<Vector3d>& points,
                               const vector<list<unsigned long>>& connectivity,
                               const vector<Matrix<double, 3, 1>>& normals,
                               double tau,
                               unsigned iter,
                               double H) {
    assert(tau > 0);
    assert(H > 0);

    auto n_pts = points.size();
    auto smoothed = points;
    for(auto i = 0; i < iter; i++) {
        auto prev_iter_pts = smoothed;
        for(auto j = 0; j < n_pts; j++) {
            auto neighbours = connectivity[j];
            auto current_point = smoothed[j];
            double total_weight = 0.0;
            Vector3d laplace = {0.0, 0.0, 0.0};
            for(auto const & neighbour : neighbours) {
                auto diff = current_point - prev_iter_pts[neighbour];
                double dst = diff.norm();
                total_weight += 1.0/dst;
                laplace += prev_iter_pts[neighbour] / dst;
            }
            smoothed[j] = prev_iter_pts[j] + tau *
                         (1.0 / total_weight * laplace - prev_iter_pts[j]);
        }

        auto max_diff = H/iter;

        for(auto j = 0; j < n_pts; j++) {
            auto diff = smoothed[j] - prev_iter_pts[j];
            auto product = diff.dot(normals[j]);
            if(product < 0) {
                product = 0.0;
                smoothed[j] = prev_iter_pts[j];
            } else {
                if(product > max_diff)
                    product = max_diff;
                smoothed[j] = prev_iter_pts[j] + product * diff / diff.norm();
            }
        }
    }
    return smoothed;
}

#include<cfenv>
int main()
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // Vis Pipeline ----------------------------
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkRenderer> renderer;

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetSize(640, 480);
    renderWindow->AddRenderer(renderer);

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(renderWindow);

    renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
    renderer->UseHiddenLineRemovalOn();
    // -------------------------------------------

    // load mesh
    std::cout << "Loading... " << std::endl;
    auto grid = ReadUnstructuredGrid("nex_nut.vtk");

//    cout << grid->GetCellType(0) << endl; // 5 is triangle
    auto n_faces = grid->GetNumberOfCells();
    auto vtkfaces = grid->GetCells();

    vector<Matrix<unsigned long, 3, 1>> faces(n_faces);

    for(auto i = 0; i < n_faces; i++) {
        auto idarr = vtkSmartPointer<vtkIdList>::New();
        vtkfaces->GetCellAtId(i, idarr);
        for(auto j = 0; j < 3; j++) // TRIANGLES ARE HARDCODED!!!!
            faces[i](j) = idarr->GetId(j);
    }
//    print(faces[0]);


    vtkSmartPointer<vtkPoints> vtkpts = grid->GetPoints();
    auto n_pts = vtkpts->GetNumberOfPoints();

    vector<Vector3d> pts(n_pts);

    for(auto i = 0; i < n_pts; i++) {
        auto p = vtkpts->GetPoint(i);
        for(auto j = 0; j < 3; j++)
            pts[i](j) = p[j];
    }

//    print(pts);

    // computations start here ----------------------------------------------
    auto normals = get_normals(pts, faces, 1);

//    // Load normals to vtk format
//    auto vtknormals = vtkSmartPointer<vtkDoubleArray>::New();
//    vtknormals->SetName("Normals");
//    vtknormals->SetNumberOfComponents(3);
//
//    for(auto i = 0; i < n_pts; i++) {
//        double arg[3] = {normals[i][0], normals[i][1], normals[i][2]};
//        vtknormals->InsertNextTuple(arg);
//    }

    // CONSTANTS
    auto m = 10;
    auto d = 0.0015;
    auto ce = 0.0008;
    double Tm = 0.0;

    auto new_layer = pts;
    auto connectivity_array = get_connectivity_array(pts.size(), faces);

    pts = vector<Vector3d>(0);

    for(auto j = 0; j < m; j++) {
        for(auto i = 0; i < n_pts; i++)
            new_layer[i] = new_layer[i] + Tm * normals[i];

        if(j != 0) {
            new_layer = laplace_points_smooth(new_layer, connectivity_array,
                                              normals, 0.5, 100, Tm);
        }

        pts.insert(pts.end(), new_layer.begin(), new_layer.end());
        Tm = pow(1 + ce, j) * d;
        normals = get_normals(new_layer, faces, 1);
        // TODO: smoothing normals
    }

    // Write to vtk the result
    auto new_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    auto new_pts = vtkSmartPointer<vtkPoints>::New();
    for(auto const & pt : pts)
        new_pts->InsertNextPoint(pt[0], pt[1], pt[2]);


    new_grid->SetPoints(new_pts);

    auto cells_size = faces.size(); // Amount of cells in 1 layer
    auto offset = normals.size(); // Amount of points in 1 layer
    for(unsigned long i = 0; i < cells_size * (m-1); i++) {
        auto layer = i / cells_size;
        auto prism = vtkSmartPointer<vtkWedge>::New();
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
        new_grid->InsertNextCell(prism->GetCellType(), prism->GetPointIds());
    }

    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName("result.vtk");
    writer->SetInputData(new_grid);
//    writer->SetFileVersion(vtkDataWriter::VTK_LEGACY_READER_VERSION_4_2);
    writer->Write();

    // Visualize ------------------------------
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(new_grid);
    mapper->ScalarVisibilityOff();

    vtkNew<vtkProperty> backProp;
    backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
    backProp->SetSpecular(.6);
    backProp->SetSpecularPower(30);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->SetBackfaceProperty(backProp);
    actor->GetProperty()->SetDiffuseColor(
            colors->GetColor3d("Tomato").GetData());
    actor->GetProperty()->SetSpecular(.3);
    actor->GetProperty()->SetSpecularPower(30);
    actor->GetProperty()->EdgeVisibilityOn();
    renderer->AddActor(actor);
    renderer->GetActiveCamera()->Azimuth(45);
    renderer->GetActiveCamera()->Elevation(45);
    renderer->ResetCamera();
    renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
    renderWindow->Render();
    interactor->Start();

    return EXIT_SUCCESS;
}
