#include "ply_utils.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>

#include "omp.h"
#include "ply.h"

namespace {

const char index_prop_name[] = "vertex_indices";

#define assert_eq(a, b) assert(a == b)
#define assert_ne(a, b) assert(a != b)
constexpr float EPS = 1e-12;
constexpr float PI = 3.14159265;

const std::vector<PlyProperty> vertex_properties{
    {const_cast<char*>("x"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, x), 0, 0, 0,
     0},
    {const_cast<char*>("y"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, y), 0, 0, 0,
     0},
    {const_cast<char*>("z"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, z), 0, 0, 0,
     0},
};

const std::vector<PlyProperty> face_properties{
    {const_cast<char*>(index_prop_name), PLY_INT, PLY_INT,
     offsetof(Face, verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(Face, size)},
};

std::vector<const char*> result_elems{"vertex", "face"};

const char result_path[] = "../src/data/result.ply";
}  // namespace

ScopedPlyFile::ScopedPlyFile(const char* path, Mode mode) {
  int type;
  float version;
  mode_ = mode;
  switch (mode) {
    case Mode::READ:
      ply_fp_ = ply_open_for_reading(const_cast<char*>(path), &n_elems_,
                                     &elems_names_, &type, &version);
      break;
    case Mode::WRITE:
      mode_ = Mode::WRITE;
      ply_fp_ =
          ply_open_for_writing(path, result_elems.size(), result_elems.data(),
                               PLY_BINARY_LE, &version);
      break;
  }
  if (!ply_fp_) {
    std::cerr << "Problems with reading file: <" << path << ">\n";
  }
}

ScopedPlyFile::~ScopedPlyFile() {
  ply_close(ply_fp_);
}

PlyModel::PlyModel(ScopedPlyFile* fp) {
  assert(!!fp);
  file_ = fp;
  PrepareModel();
}

int PlyModel::mapped_vertex(int i) const {
  while (remapping.find(i) != remapping.end())
    i = remapping.at(i);
  return i;
}

bool PlyModel::ShouldShowAfterUpdate(Face& f) {
  std::set<int> verteces;
  for (int i = 0; i < f.size; i++)
    verteces.insert(mapped_vertex(f.verts[i]));
  int new_size = verteces.size();
  if (new_size < 3) {
    f.rejected = true;
    deprecated_++;
    return false;
  }

  verteces.clear();
  auto* res = new int[new_size];
  for (int i = 0, j = 0; i < f.size; ++i) {
    const int map_i = mapped_vertex(f.verts[i]);
    if (!verteces.count(map_i)) {
      verteces.insert(map_i);
      res[j++] = map_i;
    }
  }

  f.size = new_size;
  delete[] f.verts;
  f.verts = res;
  return true;
}

void PlyModel::PrepareFaces() {
  for (auto& f : flist_)
    ShouldShowAfterUpdate(f);
}

void PlyModel::Serialize(const char* const output) {
  ScopedPlyFile file(output, ScopedPlyFile::Mode::WRITE);
  assert(!!file);

  std::vector<PlyProperty> res_props{
      {const_cast<char*>("x"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, x),
       PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      {const_cast<char*>("y"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, y),
       PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      {const_cast<char*>("z"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, z),
       PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      {const_cast<char*>("nx"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, nx),
       PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      {const_cast<char*>("ny"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, ny),
       PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      {const_cast<char*>("nz"), PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, nz),
       PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      // {const_cast<char*>("red"), PLY_UCHAR, PLY_UCHAR,
      //  offsetof(Vertex, red), PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      // {const_cast<char*>("green"), PLY_UCHAR, PLY_UCHAR,
      //  offsetof(Vertex, green), PLY_SCALAR, PLY_UINT, PLY_UINT, 0},
      // {const_cast<char*>("blue"), PLY_UCHAR, PLY_UCHAR,
      //  offsetof(Vertex, blue), PLY_SCALAR, PLY_UINT, PLY_UINT, 0}
  };

  PrepareFaces();

#pragma omp parallel for
  for (int i = 0; i < vlist_.size(); ++i) {
    auto& v = vlist_[i];
    auto& p = vertex_norm_[i].p_;
    std::tie(v.nx, v.ny, v.nz) = std::tie(p.x, p.y, p.z);
  }

  ply_describe_element(file.ply_fp(), result_elems[0], vlist_.size(),
                       res_props.size(), res_props.data());
  ply_describe_element(file.ply_fp(), result_elems[1],
                       flist_.size() - deprecated_, face_properties.size(),
                       face_properties.data());
  ply_header_complete(file.ply_fp());

  ply_put_element_setup(file.ply_fp(), result_elems[0]);

  int good = 0;
  float res = 0;
  for (auto& v : vlist_) {
    if (!v.rejected) {
      good++;
      res += v.mean_dist_to_neigbors;
    }
  }
  std::cerr << "good: " << good << " mean: " << (res / good) << '\n';

  for (auto& v : vlist_)
    ply_put_element(file.ply_fp(), (void*)(&v));

  ply_put_element_setup(file.ply_fp(), result_elems[1]);
  for (auto& f : flist_) {
    if (!f.rejected)
      ply_put_element(file.ply_fp(), (void*)(&f));
  }
}

void PlyModel::Stabilize() {
  for (int v = 0; v < vlist_.size(); v++) {
    if (vlist_[v].rejected || vlist_[v].prime)
      continue;
    float& mean = vlist_[v].mean_dist_to_neigbors;
    while (0 < mean && mean < 0.8 * density_ && collapse_edge(v)) {
    }
  }
}

void PlyModel::Simplify(float density) {
  density_ = density;
  for (int i : res_list_) {
    float& mean = vlist_[i].mean_dist_to_neigbors;
    while (0 < mean && mean < 0.8 * density_ && collapse_edge(i)) {
    }
  }
  Stabilize();
}

void PlyModel::UpdateVertexMeanDistance(int i) {
  if (!neighbor_verteces_[i].size()) {
    vlist_[i].mean_dist_to_neigbors = 0;
    return;
  }

  for (int v : neighbor_verteces_[i])
    vlist_[i].mean_dist_to_neigbors += (vlist_[v] - vlist_[i]).len();
  vlist_[i].mean_dist_to_neigbors /= neighbor_verteces_[i].size();
}

bool PlyModel::collapse_edge(int i) {
  int u = -1, v = -1;

  const int map_i = mapped_vertex(i);
  auto [val, ind] = min_cost_for_vertex(map_i);
  if (ind == -1)
    return false;

    u = map_i;
    v = mapped_vertex(ind);

  remapping[v] = u;

  neighbor_verteces_[u].erase(v);
  neighbor_verteces_[v].erase(u);

  for (int i : neighbor_verteces_[v]) {
    neighbor_verteces_[i].erase(v);
    neighbor_verteces_[i].insert(u);
    neighbor_verteces_[u].insert(i);
    UpdateVertexMeanDistance(i);
  }
  UpdateVertexMeanDistance(u);
  neighbor_verteces_[v].clear();
  UpdateVertexMeanDistance(v);
  vlist_[v].rejected = true;
  return true;
}

void PlyModel::PrepareModel() {
  for (int i = 0; i < file_->n_elems(); ++i) {
    PrepareElement(file_->elems_names()[i]);
  }

  ComputeNorms();
  DetectCornerPoints();

#pragma omp parallel for
  for (int i = 0; i < vlist_.size(); i++) {
    UpdateVertexMeanDistance(i);
  }
}

float PlyModel::cost(int a, int b) const {
  const float distance = (vlist_[a] - vlist_[b]).len();
  std::vector<int> common_faces;
  for (int i : neighbor_faces_[a]) {
    // dont account remapping :/
    if (flist_[i].includes(b))
      common_faces.push_back(i);
  }
  float curvature = 0;
  for (int i : neighbor_faces_[a]) {
    float curr_min = 1;
    for (int j : common_faces) {
      const float val = flist_[i].norm.dot_prod(flist_[j].norm);
      curr_min = std::min((1.f - val) / 2.f, curr_min);
    }
    curvature = std::max(curvature, curr_min);
  }
  return curvature * distance;
}

std::tuple<float, int> PlyModel::min_cost_for_vertex(int a) const {
  float res = 1e9;
  int ind = -1;
  for (int v : neighbor_verteces_[a]) {
    const float val = cost(a, v);
    if (val < res) {
      res = val;
      ind = v;
    }
  }

  return {res, ind};
}

void PlyModel::ComputeNorms() { 
  flist_.erase(std::remove_if(flist_.begin(), flist_.end(),
                              [](auto& f) { return f.size < 3; }),
               flist_.end());

  #pragma omp parallel for
  for (int i = 0; i < flist_.size(); ++i) {
    auto& f = flist_[i];
    Vector r{vlist_[f.verts[1]] - vlist_[f.verts[0]]},
        l{vlist_[f.verts[2]] - vlist_[f.verts[0]]};
    auto res = r.cross_prod(l);
    f.norm = res;
  }

  flist_.erase(std::remove_if(flist_.begin(), flist_.end(),
                              [](auto& f) { return f.norm.len() == 0.; }),
               flist_.end());

  for (auto& f : flist_)
    f.norm /= f.norm.len();

  for (int j = 0; j < flist_.size(); ++j) {
    for (int i = 0; i < flist_[j].size; ++i) {
      neighbor_faces_[flist_[j].verts[i]].insert(j);
      neighbor_verteces_[flist_[j].verts[(i + 1) % (flist_[j].size)]].insert(
          flist_[j].verts[i]);
      neighbor_verteces_[flist_[j].verts[(flist_[j].size + i - 1) %
                                         (flist_[j].size)]]
          .insert(flist_[j].verts[i]);
    }
  }

  #pragma omp parallel for
  for (int i = 0; i < neighbor_faces_.size(); i++) {
    if (!neighbor_faces_[i].size())
      continue;

    Vector res;
    for (const int n : neighbor_faces_[i]) {
      res += flist_[n].norm;
    }
    res /= res.len();
    vertex_norm_[i] = res;
  }
}

void PlyModel::DetectCornerPoints() {
  #pragma omp parallel for
  for (int i = 0; i < vlist_.size(); ++i) {
    if (!neighbor_verteces_[i].size())
      continue;
    std::vector<float> angles, dists;
    float max_angle = 0, min_angle = PI;
    float max_dist = 0, min_dist = 1e9;
    std::vector<int> neighbors2{i};
    for (const int n : neighbor_verteces_[i]) {
      const float angle = vertex_norm_[i].angle(vertex_norm_[n]);
      const float dist =
          Face(vlist_[i], vertex_norm_[i]).distance_to(vlist_[n]);
      angles.push_back(angle);
      dists.push_back(dist);
      max_angle = std::max(max_angle, angle);
      max_dist = std::max(max_dist, dist);
      min_angle = std::min(min_angle, angle);
      min_dist = std::min(min_dist, dist);
      neighbors2.push_back(n);
    }

    auto h_sum = [](float cumm, float cur) { return cumm + 1 / cur; };
    const float mean_angle =
        angles.size() /
        std::accumulate(angles.begin(), angles.end(), 0.f, h_sum);
    const float mean_dist =
        dists.size() / std::accumulate(dists.begin(), dists.end(), 0.f, h_sum);

    float max_angle2 = 0, min_angle2 = PI;
    float max_dist2 = 0, min_dist2 = 1e9;
    std::vector<float> angles2, dists2;
    for (const int n : neighbor_verteces_[i])
      for (const int m : neighbor_verteces_[n]) {
        if (std::find(neighbors2.begin(), neighbors2.end(), m) !=
            neighbors2.end())
          continue;
        const float angle = vertex_norm_[i].angle(vertex_norm_[m]);
        const float dist =
            Face(vlist_[i], vertex_norm_[i]).distance_to(vlist_[m]);
        angles2.push_back(angle);
        dists2.push_back(dist);
        max_angle2 = std::max(max_angle2, angle);
        max_dist2 = std::max(max_dist2, dist);
        min_angle2 = std::min(min_angle2, angle);
        min_dist2 = std::min(min_dist2, dist);
        neighbors2.push_back(m);
      }

    const float mean_angle2 =
        angles2.size() /
        std::accumulate(angles2.begin(), angles2.end(), 0.f, h_sum);
    const float mean_dist2 =
        dists2.size() /
        std::accumulate(dists2.begin(), dists2.end(), 0.f, h_sum);
    float val =
        max_dist == min_dist || max_angle == min_angle
            ? 0
            : 0.95 * (mean_dist - min_dist) / (max_dist - min_dist) +
                  0.05 * (mean_angle - min_angle) / (max_angle - min_angle);
    float val1 =
        max_dist2 == min_dist2 || max_angle2 == min_angle2
            ? 0
            : 0.95 * (mean_dist2 - min_dist2) / (max_dist2 - min_dist2) +
                  0.05 * (mean_angle2 - min_angle2) / (max_angle2 - min_angle2);
    val *= val1;
    assert(val >= 0.);
    assert(val <= 1.);
    vlist_[i].val = val;
  }

  auto lmax = [](float x, auto y) { return std::max(x, y.val); };
  const float max_val =
      std::accumulate(vlist_.begin(), vlist_.end(), 0.f, lmax);

  // add right part from vlist partion. which %? 
  for (int i = 0; i < vlist_.size(); ++i)
    if (vlist_[i].val >= 0.4 * max_val)
      res_list_.push_back(i);
}

void PlyModel::PrepareElement(char* elem_name) {
  int num_elems, nprops;
  ply_get_element_description(file_->ply_fp(), elem_name, &num_elems, &nprops);

  printf("element %s %d %d\n", elem_name, num_elems, nprops);

  if (equal_strings("vertex", elem_name)) {
    assert_eq(3, nprops);
    vlist_.resize(num_elems);
    vertex_norm_.resize(num_elems);
    neighbor_faces_.resize(num_elems);
    neighbor_verteces_.resize(num_elems);

    for (auto& p : vertex_properties)
      ply_get_property(file_->ply_fp(), elem_name, &p);

    for (auto& v : vlist_)
      ply_get_element(file_->ply_fp(), (void*)(&v));
  }

  if (equal_strings("face", elem_name)) {
    flist_.resize(num_elems);

    for (auto& p : face_properties)
      ply_get_property(file_->ply_fp(), elem_name, &p);

    for (auto& f : flist_)
      ply_get_element(file_->ply_fp(), (void*)(&f));
  }
}
