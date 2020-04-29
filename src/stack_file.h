#ifndef SEGYSTACK_STACK_FILE_H_
#define SEGYSTACK_STACK_FILE_H_

#include <map>
#include <string>

#include "mmap_file.h"
#include "stack_types.pb.h"

namespace segystack {

class StackFile {
 public:
  class SegyOptions {
   public:
    enum TraceHeaderAttribute {
      INLINE_NUMBER,
      CROSSLINE_NUMBER,
      X_COORDINATE,
      Y_COORDINATE
    };

    SegyOptions();

    void setUtmZone(int num, char zone);
    UTMZone getUtmZone() const { return utm_zone_; }

    void setTraceHeaderOffset(TraceHeaderAttribute attr, int offset);
    int getTraceHeaderOffset(TraceHeaderAttribute attr) const {
      return offsets_.at(attr);
    }

   private:
    UTMZone utm_zone_;
    std::map<TraceHeaderAttribute, int> offsets_;
  };

  StackFile(const std::string& filename);

  StackFile(const std::string& filename,
            const std::string& segy_filename,
            const SegyOptions& opts);

  const Grid& grid() const;

  int getNumInlines() const;
  int getNumCrosslines() const;

  const UTMZone& getUtmZone() const;

  void readInline(int il,
                  std::vector<float>& data,
                  float fill_value = 0.0f) const;
  void readInline(int il,
                  float* buffer,
                  size_t buffer_size,
                  float fill_value = 0.0f) const;

  void readCrossline(int xl,
                     std::vector<float>& data,
                     float fill_value = 0.0f) const;
  void readCrossline(int xl,
                     float* buffer,
                     size_t buffer_size,
                     float fill_value = 0.0f) const;

  void readDepthSlice(unsigned int sample_index,
                      std::vector<float>& data,
                      float fill_value = 0.0f) const;
  void readDepthSlice(unsigned int sample_index,
                      float* buffer,
                      size_t buffer_size,
                      float fill_value = 0.0f) const;

  bool isOptimizedForCrosslineAccess() const {
    return data_xl_file_ == nullptr;
  }
  bool isOptimizedForDepthSliceAccess() const {
    return data_ds_file_ == nullptr;
  }

  void setCrosslineAccessOptimization(bool value);
  void setDepthSliceAccessOptimization(bool value);

  ~StackFile();

 private:
  void initialize();
  void computeInlineMetadata(StackHeader::SliceMetadata* il_metadata);
  void computeCrosslineMetadata(StackHeader::SliceMetadata* xl_metadata);
  void computeDepthSliceMetadata(StackHeader::SliceMetadata* depth_metadata);
  void writeCrosslineSlices();

  void readCrosslineOptimized(int xl, float* buffer, size_t buffer_size) const;
  void readCrosslineFromInlineData(int xl,
                                   float* buffer,
                                   size_t buffer_size) const;

  void writeDepthSlices();
  void readDepthSliceOptimized(unsigned int sample_index,
                               float* buffer,
                               size_t buffer_size) const;
  void readDepthSliceFromInlineData(unsigned int sample_index,
                                    float* buffer,
                                    size_t buffer_size) const;

  class GridMap;
  std::string filename_;
  std::unique_ptr<StackHeader> header_;
  std::unique_ptr<GridMap> grid_map_;
  const Grid* grid_;
  std::unique_ptr<MmapFile> data_file_, data_xl_file_, data_ds_file_;
};

}  // namespace segystack

std::ostream& operator<<(std::ostream& os, const segystack::Grid& grid);
std::ostream& operator<<(std::ostream& os, const segystack::Grid::Cell& cell);
std::ostream& operator<<(std::ostream& os, const segystack::UTMZone& utm);
std::ostream& operator<<(std::ostream& os, const segystack::StackFile::SegyOptions& opts);

#endif  // SEGYSTACK_STACK_FILE_H_