#include "Utils.h"

Vector3d pl::get_normal_tr(const Vector3d& p1, const Vector3d& p2,
                                  const Vector3d& p3, int sign) {
    Vector3d v1 = p3-p2;
    Vector3d v2 = p1-p2;
    return v1.cross(v2) * sign;
}