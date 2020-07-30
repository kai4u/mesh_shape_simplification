#ifndef GEOM_UTILS_H_
#define GEOM_UTILS_H_

#include <iostream>

template<typename T>
T det(T a, T b, T c, T d) {
  return a * c - b * d;
}

class Vector;

struct Vertex {
  Vertex() = default;
  Vertex(float x, float y, float z) : x(x), y(y), z(z) {}
  Vector operator-(const Vertex& a) const;
  Vertex& operator*=(float a);
  Vertex& operator/=(float a);

  float x, y, z;
  float nx, ny, nz;
  unsigned char red = 0;
  unsigned char green = 0;
  unsigned char blue = 0;
  float val;
  bool prime = false;
  bool rejected = false;
  float mean_dist_to_neigbors = 0;
};

struct Vector {
  Vector() = default;
  explicit Vector(Vertex v) :p_(v) {}

  float dot_prod(const Vector& v) const;
  float angle(const Vector& v) const;
  float len() const;
  Vector cross_prod(const Vector& other) const;

  Vector& operator/=(float a);
  Vector& operator+=(const Vector& v);

  float x() const { return p_.x; }
  float y() const { return p_.y; }
  float z() const { return p_.z; }

  Vertex p_;
};

struct Face {
  Face() = default;
  Face(const Vertex& p, const Vector& norm) : p(p), norm(norm) {}

  friend std::ostream& operator<<(std::ostream& os, const Face& f) {
    os << "size: " << f.size << " elems: ";
    for (int i = 0; i < f.size; i++)
      os << f.verts[i] << ' ';
    os << '\n';
    return os;
  }

  float distance_to(const Vertex& point) const;
  bool includes(int ind) const;

  int size;
  int* verts;
  Vector norm;
  Vertex p;
  bool rejected = false;
};

#endif  // GEOM_UTILS_H_
