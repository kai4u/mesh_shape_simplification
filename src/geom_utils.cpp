#include "geom_utils.h"

#include <cassert>
#include <cmath>

namespace {

#define assert_eq(a, b) assert(a == b)
#define assert_ne(a, b) assert(a != b)

}  // namespace

Vertex& Vertex::operator*=(float a) {
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

Vertex& Vertex::operator/=(float a) {
  assert_ne(a, 0);
  return this->operator*=(1 / a);
}

Vector Vertex::operator-(const Vertex& a) const {
  return Vector(Vertex(x - a.x, y - a.y, z - a.z));
}

float Vector::dot_prod(const Vector& v) const {
  return p_.x * v.p_.x + p_.y * v.p_.y + p_.z * v.p_.z;
}

float Vector::angle(const Vector& v) const {
  float res = dot_prod(v);
  res /= len() * v.len();
  res = std::max(res, -1.f);
  res = std::min(res, 1.f);
  return acos(res);
}

float Vector::len() const {
  return sqrt(dot_prod(*this));
}

Vector Vector::cross_prod(const Vector& other) const {
  return Vector(Vertex(det(y(), z(), other.y(), other.z()),
                       det(z(), x(), other.z(), other.x()),
                       det(x(), y(), other.x(), other.y())));
}

Vector& Vector::operator/=(float a) {
  assert_ne(a, 0.);
  p_ /= a;
  return *this;
}

Vector& Vector::operator+=(const Vector& v) {
  p_.x += v.x();
  p_.y += v.y();
  p_.z += v.z();
  return *this;
}

float Face::distance_to(const Vertex& point) const {
  const float norm_len = norm.len();
  const float f_point =
      std::abs(norm.dot_prod(Vector{point}) - norm.dot_prod(Vector{p}));
  return f_point / norm_len;
}

bool Face::includes(int ind) const {
  for (int i = 0; i < size; i++)
    if (verts[i] == ind)
      return true;

  return false;
}
