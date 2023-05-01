#include "Grid.h"

using pl::Grid;

int main(int argc, char *argv[]) {
    auto grid = Grid("../data/nex_nut.vtk", 1);
    grid.constructPL(10, 0.015, 0.0008, 10, 0.5);
    grid.exportPLToVTK("result.vtk");
    grid.exportNormalsToVTK("normals.vtk");
    return 0;
}
