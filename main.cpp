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
//#include <vtkFeatureEdges.h>

// Eigen3 Headers
#include <Eigen/Dense>


using std::cout;
using std::endl;
using std::list;
using std::vector;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Vector3d;


vtkSmartPointer<vtkUnstructuredGrid> ReadUnstructuredGrid(std::string const& fileName)
{
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
void print(T& input)
{
    cout << input << endl;
}

//Matrix<double, 3, Dynamic> get_normals(
//        Matrix<double, 3, Dynamic> pc,
//        Matrix<int64_t, 3, Dynamic> faces
//        )
//{
//
//}

Vector3d get_normal_tr(const Vector3d& p1,
                       const Vector3d& p2,
                       const Vector3d& p3)
{
    auto v1 = p3-p2;
    auto v2 = p1-p2;
//    return {p1(1) * (p2(2) - p3(2)) + p2(1) * (p3(2) - p1(2)) + p3(1) * (p1(2) - p2(2)),
//            p1(0) * (p3(2) - p2(2)) + p2(0) * (p1(2) - p3(2)) + p3(0) * (p2(2) - p1(2)),
//            p1(0) * (p2(1) - p3(1))}
    return v1.cross(v2)*(-1.0);
}

int main()
{
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
    auto grid = ReadUnstructuredGrid("cube.vtk");

//    cout << grid->GetCellType(0) << endl; // 5 is triangle
    auto m = grid->GetNumberOfCells();
    auto vtkfaces = grid->GetCells();

    vector<Matrix<int64_t, 3, 1>> faces(m);

    for(auto i = 0; i < m; i++) {
        auto idarr = vtkSmartPointer<vtkIdList>::New();
        vtkfaces->GetCellAtId(i, idarr);
        for(auto j = 0; j < 3; j++) // TRIANGLES ARE HARDCODED!!!!
            faces[i](j) = idarr->GetId(j);
    }
//    print(faces[0]);


    vtkSmartPointer<vtkPoints> vtkpts = grid->GetPoints();
    auto n = vtkpts->GetNumberOfPoints();

    vector<Vector3d> pts(n);

    for(auto i = 0; i < n; i++) {
        auto p = vtkpts->GetPoint(i);
        for(auto j = 0; j < 3; j++)
            pts[i](j) = p[j];
    }

//    print(pts);

    // computations start here
    auto sub_normals = vector<list<Matrix<double, 3, 1>>>(n);
    for(auto i = 0; i < m; i++) {
        auto tri_pts = faces[i];
        auto tri_normal = get_normal_tr(pts[tri_pts(0)],
                                                 pts[tri_pts(1)],
                                                 pts[tri_pts(2)]);
        for(const auto point : faces[i]) {
//            cout << point << endl;
//            cout << tri_normal << endl;
            sub_normals[point].push_back(tri_normal);
        }
    }

    auto normals = vector<Vector3d>(n);
    auto vtknormals = vtkSmartPointer<vtkDoubleArray>::New();
    vtknormals->SetName("Normals");
    vtknormals->SetNumberOfComponents(3);

    int counter = 0;
    for(auto i = 0; i < n; i++) {
        normals[i] = {0.0, 0.0, 0.0};
        auto tmp = sub_normals[i].size();
//        print(tmp);
        counter += tmp;
        for(const auto& it : sub_normals[i]) {
//            if(it.norm() < 0.02)
//                print(it);
            normals[i] += it;
//            counter++;
        }
        if(normals[i].norm() < 0.03)
            print(normals[i]);
        normals[i] = normals[i] / normals[i].norm();
        double arg[3] = {normals[i][0], normals[i][1], normals[i][2]};
        vtknormals->InsertNextTuple(arg);
    }

    print(n);
    print(m);
    print(counter);
    print(normals[146]);
    print(normals[145]);
    print(normals[147]);
    cout << vtknormals->GetNumberOfTuples() << endl;

    grid->GetPointData()->AddArray(vtknormals);

    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName("with_normals.vtk");
    writer->SetInputData(grid);
    writer->Write();

    // Visualize ------------------------------
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(grid);
    mapper->ScalarVisibilityOff();

    vtkNew<vtkProperty> backProp;
    backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
    backProp->SetSpecular(.6);
    backProp->SetSpecularPower(30);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->SetBackfaceProperty(backProp);
    actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
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
