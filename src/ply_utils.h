#ifndef PLY_UTILS_H_
#define PLY_UTILS_H_

#include <cstdint>
#include <set>
#include <unordered_map>
#include <vector>

#include "geom_utils.h"

class PlyFile;

struct ScopedPlyFile {
  enum class Mode : uint8_t {
    READ = 0,
    WRITE,
  };

  ScopedPlyFile(const char* path, Mode mode);
  ~ScopedPlyFile();

  operator bool() const { return !!ply_fp_; }

  PlyFile* ply_fp() const { return ply_fp_; }
  int n_elems() const { return n_elems_; }
  char** elems_names() const { return elems_names_; }

 private:
  Mode mode_ = Mode::READ;
  PlyFile* ply_fp_ = nullptr;
  int n_elems_;
  char** elems_names_;
};

struct PlyModel {
  explicit PlyModel(ScopedPlyFile* fp);

  void Serialize(const char* const ouput);
  void Simplify(float density);

 private:
  bool ShouldShowAfterUpdate(Face& f);
  void PrepareFaces();
  void Stabilize();
  int mapped_vertex(int i) const;
  void UpdateVertexMeanDistance(int i);
  void PrepareModel();
  void PrepareElement(char* elem_name);

  bool collapse_edge(int i);
  float cost(int a, int b) const;
  std::tuple<float, int> min_cost_for_vertex(int a) const;

  void ComputeNorms();
  void DetectCornerPoints();


  ScopedPlyFile* file_ = nullptr;
  float density_;
  std::unordered_map<int, int> remapping;
  std::vector<Vertex> vlist_;
  std::vector<int> res_list_;
  int deprecated_ = 0; 
  std::vector<Face> flist_;
  std::vector<Vector> vertex_norm_;
  // use flat set
  std::vector<std::set<int>> neighbor_faces_;
  std::vector<std::set<int>> neighbor_verteces_;
};

#endif  // PLY_UTILS_H_
