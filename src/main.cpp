#include <chrono>
#include <cstdlib>
#include <iostream>

#include "geom_utils.h"
#include "ply.h"
#include "ply_utils.h"

namespace {

struct ScopedMeasure {
  ScopedMeasure() { start_ = std::chrono::steady_clock::now(); }
  ~ScopedMeasure() {
    const auto stop = std::chrono::steady_clock::now();
    const auto dur = std::chrono::duration_cast<std::chrono::duration<double>>(
        stop - start_);
    std::cerr << "Time: " << dur.count() << " seconds\n";
  }

 private:
  std::chrono::steady_clock::time_point start_;
};

// const char model_path[] = "../src/data/citywall_3kk.ply";
// const char model_path[] = "../src/data/unicorn0.ply";

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "call ./main <input_file> <output_file> density\n";
    return 1;
  }
  if (ScopedPlyFile file{argv[1], ScopedPlyFile::Mode::READ}) {
    ScopedMeasure measure;
    PlyModel model(&file);
    model.Simplify(atof(argv[3]));
    model.Serialize(argv[2]);
  }
  return 0;
}
