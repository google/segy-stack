/*
 * Copyright 2020 Google LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SEGYSTACK_STACK_FILE_H_
#define SEGYSTACK_STACK_FILE_H_

#include <map>
#include <string>
#include <utility>

#include "crs_util.h"
#include "mmap_file.h"
#include "segy_file.h"
#include "stack_types.pb.h"

namespace segystack {

class StackFile {
  class GridMap;

 public:
  class Grid {
   public:
    enum Units { METERS = 0, FEET = 1 };

    Grid();
    ~Grid();

    uint32_t numInlines() const;
    void setNumInlines(uint32_t value);

    uint32_t numCrosslines() const;
    void setNumCrosslines(uint32_t value);

    UTMZone utmZone() const;

    Units units() const;
    void setUnits(Units value);

    // Smallest inline number.
    int32_t inlineMin() const;
    void setInlineMin(int32_t value);

    // Largest inline number.
    int32_t inlineMax() const;
    void setInlineMax(int32_t value);

    // Inline increment.
    uint32_t inlineIncrement() const;
    void setInlineIncrement(uint32_t value);

    // Inline spacing.
    float inlineSpacing() const;
    void setInlineSpacing(float value);

    // Smallest crossline number.
    int32_t crosslineMin() const;
    void setCrosslineMin(int32_t value);

    // Largest crossline number.
    int32_t crosslineMax() const;
    void setCrosslineMax(int32_t value);

    // Crossline increment.
    uint32_t crosslineIncrement() const;
    void setCrosslineIncrement(uint32_t value);

    // Crossline spacing.
    float crosslineSpacing() const;
    void setCrosslineSpacing(float value);

    // Depth sampling interval.
    float samplingInterval() const;
    void setSamplingInterval(float value);

    // Number of samples in depth.
    uint32_t numSamples() const;
    void setNumSamples(uint32_t value);

    bool isBinActive(int32_t inline_num, int32_t crossline_num) const;

    struct Coordinate {
      float x, y;
      int32_t inline_num, crossline_num;
      float lat, lon;
      Coordinate()
          : x(0), y(0), inline_num(0), crossline_num(0), lat(0), lon(0) {}
      Coordinate(float x, float y)
          : x(x), y(y), inline_num(0), crossline_num(0), lat(0), lon(0) {}
      Coordinate(float x, float y, int32_t il, int32_t xl)
          : x(x), y(y), inline_num(il), crossline_num(xl), lat(0), lon(0) {}
      Coordinate(float x, float y, int32_t il, int32_t xl, float lat, float lon)
          : x(x), y(y), inline_num(il), crossline_num(xl), lat(lat), lon(lon) {}
    };

    Coordinate getCoordinate(int32_t inline_num, int32_t crossline_num) const;

    struct BoundingBox {
      Coordinate c1, c2, c3, c4;
    };

    BoundingBox boundingBox() const;

   protected:
    friend class StackFile;
    friend std::ostream& operator<<(std::ostream& os,
                                    const segystack::StackFile::Grid& grid);

    Grid(const GridMap* grid_map, const UTMZone& utm);
    const GridMap* grid_map_;
    UTMZone utm_zone_;
    internal::GridData grid_data_;
    std::unique_ptr<UTMZoneConverter> utm_converter_;
  };

  class SegyOptions {
   public:
    SegyOptions();

    void setUtmZone(int num, char zone);
    UTMZone getUtmZone() const { return utm_zone_; }

    void setTraceHeaderOffset(SegyFile::Trace::Header::Attribute attr,
                              int offset);
    int getTraceHeaderOffset(SegyFile::Trace::Header::Attribute attr) const {
      return offsets_.at(attr);
    }

    void setTraceHeaderOffsets(
        const std::map<SegyFile::Trace::Header::Attribute, int>& offsets);

    bool is2D() const { return is_2d_; }
    void setIs2D(bool value) { is_2d_ = value; }

   private:
    UTMZone utm_zone_;
    std::map<SegyFile::Trace::Header::Attribute, int> offsets_;
    bool is_2d_;
  };

  explicit StackFile(const std::string& filename);

  explicit StackFile(const std::string& filename, const SegyFile& segyfile);

  explicit StackFile(const std::string& filename,
                     const SegyFile& segyfile,
                     const SegyOptions& opts);

  explicit StackFile(const SegyFile& segyfile, const SegyOptions& opts);

  const std::string& name() const { return filename_; }

  const Grid& grid() const;

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
  void createFromSegy(const std::string& filename,
                      const SegyFile& segyfile,
                      const SegyOptions& opts);
  void initialize();
  void computeInlineMetadata(internal::StackHeader::SliceMetadata* il_metadata);
  void computeCrosslineMetadata(
      internal::StackHeader::SliceMetadata* xl_metadata);
  void computeDepthSliceMetadata(
      internal::StackHeader::SliceMetadata* depth_metadata);
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

  std::string filename_;
  std::unique_ptr<internal::StackHeader> header_;
  std::unique_ptr<GridMap> grid_map_;
  std::unique_ptr<const Grid> grid_;
  std::unique_ptr<MmapFile> data_file_, data_xl_file_, data_ds_file_;
};

std::ostream& operator<<(std::ostream& os, const StackFile::Grid& grid);
std::ostream& operator<<(std::ostream& os,
                         const StackFile::Grid::Coordinate& coord);
std::ostream& operator<<(std::ostream& os,
                         const StackFile::Grid::BoundingBox& bbox);
std::ostream& operator<<(std::ostream& os,
                         const StackFile::Grid::Coordinate& coord);
std::ostream& operator<<(std::ostream& os, const StackFile::SegyOptions& opts);

}  // namespace segystack

#endif  // SEGYSTACK_STACK_FILE_H_