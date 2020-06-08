// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "stack_file.h"

#include <glog/logging.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <unordered_map>

#include "logging.h"
#include "segy_file.h"

using segystack::internal::GridData;

namespace segystack {

std::ostream& operator<<(std::ostream& os, const GridData& grid_data) {
  os << "IL range : " << grid_data.inline_min() << " - "
     << grid_data.inline_max() << std::endl;
  os << "XL range : " << grid_data.crossline_min() << " - "
     << grid_data.crossline_max() << std::endl;
  os << "IL increment : " << grid_data.inline_increment() << std::endl;
  os << "XL increment : " << grid_data.crossline_increment() << std::endl;
  os << "IL spacing : " << grid_data.inline_spacing() << std::endl;
  os << "XL spacing : " << grid_data.crossline_spacing() << std::endl;
  os << "Sampling interval : " << grid_data.sampling_interval() << std::endl;
  os << "Num depth samples : " << grid_data.num_samples() << std::endl;
  os << "Units : ";
  switch (grid_data.units()) {
    case GridData::METERS:
      os << "Meters";
      break;
    case GridData::FEET:
      os << "Feet";
      break;
    default:
      os << "Unknown";
  };
  os << std::endl;
  os << "Num active cells : " << grid_data.num_active_cells() << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const internal::UTMZone& utm) {
  os << utm.number() << " " << utm.letter();
  return os;
}

std::ostream& operator<<(std::ostream& os, const StackFile::Grid& grid) {
  os << "IL num: " << grid.numInlines() << std::endl;
  os << "XL num: " << grid.numCrosslines() << std::endl;
  os << grid.grid_data_;
  os << "UTM Zone: " << grid.utm_zone_ << std::endl;
  os << "Bounding box: " << grid.boundingBox() << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const StackFile::Grid::Coordinate& coord) {
  os << "X/Y(" << coord.x << ", " << coord.y << "), IL/XL(" << coord.inline_num
     << ", " << coord.crossline_num << "), Lat/Lon(" << coord.lat << ", "
     << coord.lon << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const StackFile::Grid::BoundingBox& bbox) {
  os << std::endl << "C3: " << bbox.c3 << "\t---\tC4: " << bbox.c4;
  os << std::endl << "C1: " << bbox.c1 << "\t---\tC2: " << bbox.c2 << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const GridData::Cell& cell) {
  os << "Coordinate (x, y) : (" << cell.x_coordinate() << ", "
     << cell.y_coordinate() << ")" << std::endl;
  os << "Grid numbers (IL, XL) : (" << cell.inline_number() << ", "
     << cell.crossline_number() << ")" << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const StackFile::SegyOptions& opts) {
  using Attribute = SegyFile::Trace::Header::Attribute;
  os << "UTM zone: " << opts.getUtmZone() << std::endl;
  os << "Inline number offset: "
     << opts.getTraceHeaderOffset(Attribute::INLINE_NUMBER) << std::endl;
  os << "Crossline number offset: "
     << opts.getTraceHeaderOffset(Attribute::CROSSLINE_NUMBER) << std::endl;
  os << "X Coordinate offset: "
     << opts.getTraceHeaderOffset(Attribute::X_COORDINATE) << std::endl;
  os << "Y Coordinate offset: "
     << opts.getTraceHeaderOffset(Attribute::Y_COORDINATE) << std::endl;
  os << "Is 2D: " << opts.is2D() << std::endl;
  return os;
}

constexpr char kStackFileMagic[] = "**StackFile**\n";
constexpr int kStackFileVersion = 1;
constexpr int kSegyOffsetSampleInterval = 17;
constexpr int kSegyOffsetNumSamples = 21;
constexpr int kSegyOffsetMeasurementUnits = 55;

StackFile::~StackFile() = default;

StackFile::SegyOptions::SegyOptions() : is_2d_(false) {
  // Set defaults based on the SEGY standard.
  offsets_[SegyFile::Trace::Header::Attribute::X_COORDINATE] = 181;
  offsets_[SegyFile::Trace::Header::Attribute::Y_COORDINATE] = 185;
  offsets_[SegyFile::Trace::Header::Attribute::INLINE_NUMBER] = 189;
  offsets_[SegyFile::Trace::Header::Attribute::CROSSLINE_NUMBER] = 193;
}

void StackFile::SegyOptions::setTraceHeaderOffsets(
    const std::map<SegyFile::Trace::Header::Attribute, int>& offsets) {
  for (const auto& it : offsets) {
    setTraceHeaderOffset(it.first, it.second);
  }
}

void StackFile::SegyOptions::setUtmZone(int num, char zone) {
  utm_zone_.setValue(num, zone);
}

void StackFile::SegyOptions::setTraceHeaderOffset(
    SegyFile::Trace::Header::Attribute attr,
    int offset) {
  if (offset < 1 || offset > 231)
    throw std::runtime_error("Offset specified outside of range [1, 231]");
  offsets_[attr] = offset;
}

class StackFile::GridMap {
 public:
  GridMap() {}

  GridMap(const GridData& grid_data) : grid_data_(grid_data) {}

  GridData::Cell* addCell() {
    cells_.push_back(GridData::Cell());
    return &(cells_.back());
  }

  void RecomputeGrid() {
    grid_data_.set_inline_min(std::numeric_limits<int>::max());
    grid_data_.set_crossline_min(std::numeric_limits<int>::max());
    grid_data_.set_inline_max(std::numeric_limits<int>::min());
    grid_data_.set_crossline_max(std::numeric_limits<int>::min());

    for (const GridData::Cell& cell : cells_) {
      grid_data_.set_inline_min(
          std::min(grid_data_.inline_min(), cell.inline_number()));
      grid_data_.set_crossline_min(
          std::min(grid_data_.crossline_min(), cell.crossline_number()));
      grid_data_.set_inline_max(
          std::max(grid_data_.inline_max(), cell.inline_number()));
      grid_data_.set_crossline_max(
          std::max(grid_data_.crossline_max(), cell.crossline_number()));
    }
    grid_data_.set_num_active_cells(cells_.size());

    grid_data_.set_inline_increment(computeIncrement(
        grid_data_.inline_min(), grid_data_.inline_max(),
        [](const GridData::Cell& cell) { return cell.inline_number(); }));

    grid_data_.set_crossline_increment(computeIncrement(
        grid_data_.crossline_min(), grid_data_.crossline_max(),
        [](const GridData::Cell& cell) { return cell.crossline_number(); }));

    Finalize();

    grid_data_.set_inline_spacing(
        computeSpacing(&GridMap::getCoordInNextInline));
    grid_data_.set_crossline_spacing(
        computeSpacing(&GridMap::getCoordInNextCrossline));
  }

  void Finalize() {
    cell_map_.resize(getNumInlines());
    for (int i = 0; i < getNumInlines(); i++) {
      cell_map_[i].resize(getNumCrosslines(), nullptr);
    }

    for (const GridData::Cell& cell : cells_) {
      int il_idx = getInlineIdx(cell.inline_number());
      int xl_idx = getCrosslineIdx(cell.crossline_number());
      cell_map_[il_idx][xl_idx] = &cell;
    }

    active_il_index_map_.clear();
    int il_idx = 0;
    for (int il = grid_data_.inline_min(); il <= grid_data_.inline_max();
         il += grid_data_.inline_increment()) {
      for (int xl = grid_data_.crossline_min();
           xl <= grid_data_.crossline_max();
           xl += grid_data_.crossline_increment()) {
        if (isCellActive(il, xl)) {
          active_il_index_map_[il] = il_idx++;
          break;
        }
      }
    }

    active_xl_index_map_.clear();
    int xl_idx = 0;
    for (int xl = grid_data_.crossline_min(); xl <= grid_data_.crossline_max();
         xl += grid_data_.crossline_increment()) {
      for (int il = grid_data_.inline_min(); il <= grid_data_.inline_max();
           il += grid_data_.inline_increment()) {
        if (isCellActive(il, xl)) {
          active_xl_index_map_[xl] = xl_idx++;
          break;
        }
      }
    }
  }

  int getNumInlines() const {
    CHECK_GT(grid_data_.inline_increment(), 0);
    return ((grid_data_.inline_max() - grid_data_.inline_min()) /
            grid_data_.inline_increment()) +
           1;
  }

  int getNumCrosslines() const {
    CHECK_GT(grid_data_.crossline_increment(), 0);
    return ((grid_data_.crossline_max() - grid_data_.crossline_min()) /
            grid_data_.crossline_increment()) +
           1;
  }

  int getInlineIncrement() const { return grid_data_.inline_increment(); }
  int getCrosslineIncrement() const { return grid_data_.crossline_increment(); }

  const GridData& gridData() const { return grid_data_; }

  const std::vector<GridData::Cell>& cells() const { return cells_; }

  bool isCellActive(int il, int xl) const {
    try {
      return cell_map_.at(getInlineIdx(il)).at(getCrosslineIdx(xl)) != nullptr;
    } catch (const std::out_of_range& e) {
      DLOG(INFO) << __FUNCTION__ << ": (" << il << ", " << xl
                 << ") out of bounds";
    }
    return false;
  }

  const GridData::Cell* getCell(int il, int xl) const {
    try {
      return cell_map_.at(getInlineIdx(il)).at(getCrosslineIdx(xl));
    } catch (const std::out_of_range& e) {
      DLOG(INFO) << __FUNCTION__ << ": (" << il << ", " << xl
                 << ") out of bounds";
    }
    return nullptr;
  }

  void setInlineDataAddress(const char* trc_addr) {
    size_t num_trace_bytes = sizeof(float) * grid_data_.num_samples();

    trace_map_.resize(getNumInlines());
    for (size_t i = 0; i < trace_map_.size(); i++) {
      trace_map_[i].resize(getNumCrosslines(), nullptr);
    }

    for (int il = grid_data_.inline_min(); il <= grid_data_.inline_max();
         il += grid_data_.inline_increment()) {
      for (int xl = grid_data_.crossline_min();
           xl <= grid_data_.crossline_max();
           xl += grid_data_.crossline_increment()) {
        if (isCellActive(il, xl)) {
          trace_map_[getInlineIdx(il)][getCrosslineIdx(xl)] =
              reinterpret_cast<const float*>(trc_addr);
          trc_addr += num_trace_bytes;
        }
      }
    }
  }

  const float* getTrace(int il, int xl) const {
    try {
      return trace_map_.at(getInlineIdx(il)).at(getCrosslineIdx(xl));
    } catch (const std::out_of_range& e) {
      DLOG(INFO) << __FUNCTION__ << ": (" << il << ", " << xl
                 << ") out of bounds";
    }
    return nullptr;
  }

  int getInlineIdx(int il) const {
    return (il - grid_data_.inline_min()) / grid_data_.inline_increment();
  }

  int getActiveInlineIdx(int il) const {
    if (active_il_index_map_.find(il) != active_il_index_map_.end())
      return active_il_index_map_.at(il);
    return -1;
  }

  int getCrosslineIdx(int xl) const {
    return (xl - grid_data_.crossline_min()) / grid_data_.crossline_increment();
  }

  int getActiveCrosslineIdx(int xl) const {
    if (active_xl_index_map_.find(xl) != active_xl_index_map_.end())
      return active_xl_index_map_.at(xl);
    return -1;
  }

  Grid::BoundingBox computeBoundingBox() const {
    Grid::BoundingBox bbox;

    LOGFN_VAR(getNumInlines());
    LOGFN_VAR(getNumCrosslines());

    if (getNumInlines() == 1 && getNumCrosslines() == 1) {
      // Single cell case.
      const GridData::Cell* cell = cell_map_[0][0];
      CHECK_NOTNULL(cell);
      Grid::Coordinate c(cell->x_coordinate(), cell->y_coordinate(),
                         grid_data_.inline_min(), grid_data_.crossline_min());
      bbox.c1 = c;
      bbox.c2 = c;
      bbox.c3 = c;
      bbox.c4 = c;
      return bbox;
    } else if (getNumInlines() == 1) {
      // Single inline case.
      const GridData::Cell* cell1 = cell_map_[0][0];
      const GridData::Cell* cell2 = cell_map_[0][getNumCrosslines() - 1];
      CHECK_NOTNULL(cell1);
      CHECK_NOTNULL(cell2);
      Grid::Coordinate c1(cell1->x_coordinate(), cell1->y_coordinate(),
                          grid_data_.inline_min(), grid_data_.crossline_min());
      Grid::Coordinate c2(cell2->x_coordinate(), cell2->y_coordinate(),
                          grid_data_.inline_min(), grid_data_.crossline_max());
      bbox.c1 = c1;
      bbox.c3 = c1;
      bbox.c2 = c2;
      bbox.c4 = c2;
      return bbox;
    } else if (getNumCrosslines() == 1) {
      // Single crossline case.
      const GridData::Cell* cell1 = cell_map_[0][0];
      const GridData::Cell* cell2 = cell_map_[getNumInlines() - 1][0];
      CHECK_NOTNULL(cell1);
      CHECK_NOTNULL(cell2);
      Grid::Coordinate c1(cell1->x_coordinate(), cell1->y_coordinate(),
                          grid_data_.inline_min(), grid_data_.crossline_min());
      Grid::Coordinate c3(cell2->x_coordinate(), cell2->y_coordinate(),
                          grid_data_.inline_max(), grid_data_.crossline_min());
      bbox.c1 = c1;
      bbox.c3 = c3;
      bbox.c2 = c1;
      bbox.c4 = c3;
      return bbox;
    }

    float max_len_il = 0;
    size_t origin_il_idx = 0, origin_xl_idx = 0;
    const GridData::Cell* origin = nullptr;
    const GridData::Cell* support = nullptr;
    for (size_t i = 0; i < cell_map_.size(); i++) {
      const GridData::Cell* c1 = nullptr;
      const GridData::Cell* c2 = nullptr;
      size_t c1_xl_idx = 0;
      for (size_t j = 0; j < cell_map_[i].size(); j++) {
        if (cell_map_[i][j]) {
          c1 = cell_map_[i][j];
          c1_xl_idx = j;
          break;
        }
      }

      for (int j = cell_map_[i].size() - 1; j >= 0; j--) {
        if (cell_map_[i][j]) {
          c2 = cell_map_[i][j];
          break;
        }
      }

      if (c1 && c2) {
        float dist = computeDistBetweenCells(*c1, *c2);
        if (dist > max_len_il) {
          max_len_il = dist;
          origin_il_idx = i;
          origin_xl_idx = c1_xl_idx;
          origin = c1;
          support = c2;
        }
      }
    }

    Grid::Coordinate c1_r, c2_r, c3_r, c4_r;
    float theta = 0.0f;

    if (origin == nullptr || support == nullptr ||
        computeDistBetweenCells(*origin, *support) == 0.0f) {
      LOG(WARNING) << "Warning: The traces are missing coordinate information. "
                      "The bounding box will not have coordinate information."
                   << std::endl;
    } else {
      CHECK_NOTNULL(origin);
      CHECK_NOTNULL(support);

      LOGFN_VAR(*origin);
      LOGFN_VAR(*support);

      float x_dist = support->x_coordinate() - origin->x_coordinate();
      float y_dist = support->y_coordinate() - origin->y_coordinate();

      CHECK_NE(std::abs(x_dist) == 0.0f && std::abs(y_dist) == 0.0f, true);

      theta = std::atan2(y_dist, x_dist);

      // The four corners in the local coordinate system around the origin.
      c1_r.x = -1.0 * origin_xl_idx * grid_data_.crossline_spacing();
      c1_r.y = -1.0 * origin_il_idx * grid_data_.inline_spacing();

      c2_r.x = (getNumCrosslines() - origin_xl_idx - 1) *
               grid_data_.crossline_spacing();
      c2_r.y = c1_r.y;

      c3_r.x = c1_r.x;
      c3_r.y =
          (getNumInlines() - origin_il_idx - 1) * grid_data_.inline_spacing();

      c4_r.x = c2_r.x;
      c4_r.y = c3_r.y;
    }

    auto set_corner_value = [&](int32_t il, int32_t xl,
                                const Grid::Coordinate& c_r,
                                const GridData::Cell* origin, float theta,
                                Grid::Coordinate* corner) {
      const GridData::Cell* cell = getCell(il, xl);
      if (cell) {
        corner->x = cell->x_coordinate();
        corner->y = cell->y_coordinate();
      } else if (origin) {
        (*corner) = ConvertToGlobalCRS(c_r, *origin, theta);
      }
      corner->inline_num = il;
      corner->crossline_num = xl;
    };

    set_corner_value(grid_data_.inline_min(), grid_data_.crossline_min(), c1_r,
                     origin, theta, &bbox.c1);
    set_corner_value(grid_data_.inline_min(), grid_data_.crossline_max(), c2_r,
                     origin, theta, &bbox.c2);
    set_corner_value(grid_data_.inline_max(), grid_data_.crossline_min(), c3_r,
                     origin, theta, &bbox.c3);
    set_corner_value(grid_data_.inline_max(), grid_data_.crossline_max(), c4_r,
                     origin, theta, &bbox.c4);

    return bbox;
  }

 private:
  static Grid::Coordinate ConvertToGlobalCRS(const Grid::Coordinate& c_r,
                                             const GridData::Cell& origin,
                                             float theta) {
    Grid::Coordinate c(c_r);
    // Rotate in opposite direction and then translate back.
    c.x = c_r.x * std::cos(-theta) + c_r.y * std::sin(-theta) +
          origin.x_coordinate();
    c.y = -c_r.x * std::sin(-theta) + c_r.y * std::cos(-theta) +
          origin.y_coordinate();
    return c;
  }

  int computeIncrement(int min_value,
                       int max_value,
                       std::function<int(const GridData::Cell&)> line_number) {
    if (min_value == max_value) {
      return 1;
    }

    std::vector<bool> line_number_map((max_value - min_value) + 1, false);
    for (size_t i = 0; i < cells_.size(); i++) {
      line_number_map[line_number(cells_[i]) - min_value] = true;
    }

    int smallest_increment = std::numeric_limits<int>::max();
    int prev_active_idx = -1;
    for (size_t i = 0; i < line_number_map.size(); i++) {
      if (line_number_map[i]) {
        if (prev_active_idx >= 0) {
          smallest_increment =
              std::min(smallest_increment, int(i) - prev_active_idx);
        }
        prev_active_idx = i;
      }
    }

    if (smallest_increment == std::numeric_limits<int>::max()) {
      return 1;
    }
    return smallest_increment;
  }

  typedef const GridData::Cell* (GridMap::*GetNextCoordMethod)(size_t,
                                                               size_t) const;
  float computeSpacing(GetNextCoordMethod get_next_coord) const {
    std::vector<float> spacings;
    for (size_t i = 0; i < cell_map_.size(); i++) {
      for (size_t j = 0; j < cell_map_[i].size(); j++) {
        const GridData::Cell* c1 = cell_map_[i][j];
        const GridData::Cell* c2 = (this->*get_next_coord)(i, j);
        if (c1 && c2) {
          float dist = computeDistBetweenCells(*c1, *c2);
          spacings.push_back(dist);
        }
      }
    }

    if (spacings.size() == 0) {
      return 1.0;
    }

    size_t center_idx = spacings.size() / 2;
    std::nth_element(spacings.begin(), spacings.begin() + center_idx,
                     spacings.end());
    return spacings[center_idx];
  }

  float computeDistBetweenCells(const GridData::Cell& c1,
                                const GridData::Cell& c2) const {
    return std::sqrt(std::pow(c1.x_coordinate() - c2.x_coordinate(), 2.0) +
                     std::pow(c1.y_coordinate() - c2.y_coordinate(), 2.0));
  }

  const GridData::Cell* getCoordInNextInline(size_t il_idx,
                                             size_t xl_idx) const {
    if (il_idx >= 0 && il_idx < cell_map_.size() - 1 && xl_idx >= 0 &&
        xl_idx < cell_map_[il_idx + 1].size()) {
      return cell_map_[il_idx + 1][xl_idx];
    }
    return nullptr;
  }

  const GridData::Cell* getCoordInNextCrossline(size_t il_idx,
                                                size_t xl_idx) const {
    if (il_idx >= 0 && il_idx < cell_map_.size() && xl_idx >= 0 &&
        xl_idx < cell_map_[il_idx].size() - 1) {
      return cell_map_[il_idx][xl_idx + 1];
    }
    return nullptr;
  }

  GridData grid_data_;
  std::vector<GridData::Cell> cells_;
  std::vector<std::vector<const GridData::Cell*>> cell_map_;
  std::unordered_map<int, int> active_il_index_map_, active_xl_index_map_;
  std::vector<std::vector<const float*>> trace_map_;
};

StackFile::Grid::Grid() {
  grid_data_.set_inline_min(0);
  grid_data_.set_inline_max(0);
  grid_data_.set_inline_increment(1);
  grid_data_.set_crossline_min(0);
  grid_data_.set_crossline_max(0);
  grid_data_.set_crossline_increment(1);
  utm_converter_.reset(new UTMZoneConverter(utm_zone_));
}

StackFile::Grid::Grid(const GridMap* grid_map, const UTMZone& utm)
    : grid_map_(grid_map),
      utm_zone_(utm),
      grid_data_(grid_map->gridData()),
      utm_converter_(new UTMZoneConverter(utm)) {}

StackFile::Grid::~Grid() = default;

uint32_t StackFile::Grid::numInlines() const {
  CHECK_GT(grid_data_.inline_increment(), 0);
  return ((grid_data_.inline_max() - grid_data_.inline_min()) /
          grid_data_.inline_increment()) +
         1;
}

void StackFile::Grid::setNumInlines(uint32_t value) {
  grid_data_.set_inline_max(grid_data_.inline_min() +
                            (value - 1) * grid_data_.inline_increment());
}

uint32_t StackFile::Grid::numCrosslines() const {
  CHECK_GT(grid_data_.crossline_increment(), 0);
  return std::abs((grid_data_.crossline_max() - grid_data_.crossline_min()) /
                  grid_data_.crossline_increment()) +
         1;
}

void StackFile::Grid::setNumCrosslines(uint32_t value) {
  grid_data_.set_crossline_max(grid_data_.crossline_min() +
                               (value - 1) * grid_data_.crossline_increment());
}

UTMZone StackFile::Grid::utmZone() const {
  return utm_zone_;
}

StackFile::Grid::Units StackFile::Grid::units() const {
  return static_cast<Units>(grid_data_.units());
}

void StackFile::Grid::setUnits(StackFile::Grid::Units value) {
  grid_data_.set_units(static_cast<internal::GridData_Units>(value));
}

int32_t StackFile::Grid::inlineMin() const {
  return grid_data_.inline_min();
}

void StackFile::Grid::setInlineMin(int32_t value) {
  grid_data_.set_inline_min(value);
}

int32_t StackFile::Grid::inlineMax() const {
  return grid_data_.inline_max();
}

void StackFile::Grid::setInlineMax(int32_t value) {
  CHECK_GE(value, grid_data_.inline_min());
  grid_data_.set_inline_max(value);
}

uint32_t StackFile::Grid::inlineIncrement() const {
  return grid_data_.inline_increment();
}

void StackFile::Grid::setInlineIncrement(uint32_t value) {
  CHECK_GT(value, 0u);
  grid_data_.set_inline_increment(value);
}

float StackFile::Grid::inlineSpacing() const {
  return grid_data_.inline_spacing();
}

void StackFile::Grid::setInlineSpacing(float value) {
  CHECK_GT(value, 0.0f);
  grid_data_.set_inline_spacing(value);
}

int32_t StackFile::Grid::crosslineMin() const {
  return grid_data_.crossline_min();
}

void StackFile::Grid::setCrosslineMin(int32_t value) {
  grid_data_.set_crossline_min(value);
}

int32_t StackFile::Grid::crosslineMax() const {
  return grid_data_.crossline_max();
}

void StackFile::Grid::setCrosslineMax(int32_t value) {
  CHECK_GE(value, grid_data_.crossline_min());
  grid_data_.set_crossline_max(value);
}

uint32_t StackFile::Grid::crosslineIncrement() const {
  return grid_data_.crossline_increment();
}

void StackFile::Grid::setCrosslineIncrement(uint32_t value) {
  CHECK_GT(value, 0u);
  grid_data_.set_crossline_increment(value);
}

float StackFile::Grid::crosslineSpacing() const {
  return grid_data_.crossline_spacing();
}

void StackFile::Grid::setCrosslineSpacing(float value) {
  CHECK_GT(value, 0.0f);
  grid_data_.set_crossline_spacing(value);
}

float StackFile::Grid::samplingInterval() const {
  return grid_data_.sampling_interval();
}

void StackFile::Grid::setSamplingInterval(float value) {
  CHECK_GT(value, 0.0f);
  grid_data_.set_sampling_interval(value);
}

uint32_t StackFile::Grid::numSamples() const {
  return grid_data_.num_samples();
}

void StackFile::Grid::setNumSamples(uint32_t value) {
  CHECK_GT(value, 0u);
  grid_data_.set_num_samples(value);
}

bool StackFile::Grid::isBinActive(int32_t inline_num,
                                  int32_t crossline_num) const {
  if (grid_map_)
    return grid_map_->isCellActive(inline_num, crossline_num);

  return false;
}

StackFile::Grid::Coordinate StackFile::Grid::getCoordinate(
    int32_t inline_num,
    int32_t crossline_num) const {
  if (grid_map_) {
    const GridData::Cell* cell = grid_map_->getCell(inline_num, crossline_num);
    if (cell) {
      GeographicCoordinates geo_coord =
          utm_converter_->getGeographicCoordinates(cell->x_coordinate(),
                                                   cell->y_coordinate());

      Coordinate coord(cell->x_coordinate(), cell->y_coordinate(), inline_num,
                       crossline_num, geo_coord.latitude, geo_coord.longitude);
      return coord;
    }
  }

  return Coordinate();
}

StackFile::Grid::BoundingBox StackFile::Grid::boundingBox() const {
  if (grid_map_) {
    BoundingBox bbox = grid_map_->computeBoundingBox();
    auto assign_lat_lon = [&](Coordinate* coord) {
      if (!utm_converter_)
        return;
      GeographicCoordinates geo_coord =
          utm_converter_->getGeographicCoordinates(coord->x, coord->y);
      coord->lat = geo_coord.latitude;
      coord->lon = geo_coord.longitude;
    };
    assign_lat_lon(&bbox.c1);
    assign_lat_lon(&bbox.c2);
    assign_lat_lon(&bbox.c3);
    assign_lat_lon(&bbox.c4);
    return bbox;
  }

  return BoundingBox();
}

void StackFile::computeInlineMetadata(
    internal::StackHeader::SliceMetadata* il_metadata) {
  CHECK_NOTNULL(grid_map_);

  size_t num_bytes_written = 0;
  const GridData& gd = grid_map_->gridData();

  for (int il = gd.inline_min(); il <= gd.inline_max();
       il += gd.inline_increment()) {
    size_t il_size = 0;

    for (int xl = gd.crossline_min(); xl <= gd.crossline_max();
         xl += gd.crossline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        il_size += sizeof(float) * gd.num_samples();
      }
    }

    if (il_size > 0) {
      il_metadata->add_size(il_size);
      il_metadata->add_offset(num_bytes_written);
    }
    num_bytes_written += il_size;
  }
}

void StackFile::computeCrosslineMetadata(
    internal::StackHeader::SliceMetadata* xl_metadata) {
  CHECK_NOTNULL(grid_map_);

  size_t num_bytes_written = 0;
  const GridData& gd = grid_map_->gridData();

  for (int xl = gd.crossline_min(); xl <= gd.crossline_max();
       xl += gd.crossline_increment()) {
    size_t xl_size = 0;

    for (int il = gd.inline_min(); il <= gd.inline_max();
         il += gd.inline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        xl_size += sizeof(float) * gd.num_samples();
      }
    }

    if (xl_size > 0) {
      xl_metadata->add_size(xl_size);
      xl_metadata->add_offset(num_bytes_written);
    }

    num_bytes_written += xl_size;
  }
}

void StackFile::computeDepthSliceMetadata(
    internal::StackHeader::SliceMetadata* depth_metadata) {
  CHECK_NOTNULL(grid_map_);

  size_t depth_slice_bytes = 0;
  const GridData& gd = grid_map_->gridData();

  for (int il = gd.inline_min(); il <= gd.inline_max();
       il += gd.inline_increment()) {
    for (int xl = gd.crossline_min(); xl <= gd.crossline_max();
         xl += gd.crossline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        depth_slice_bytes += sizeof(float);
      }
    }
  }

  size_t num_bytes_written = 0;
  for (size_t iz = 0; iz < gd.num_samples(); iz++) {
    depth_metadata->add_offset(num_bytes_written);
    num_bytes_written += depth_slice_bytes;
    depth_metadata->add_size(depth_slice_bytes);
  }
}

StackFile::StackFile(const std::string& filename) : filename_(filename) {
  {
    MmapFile fp(filename_);
    if (!fp.exists()) {
      throw std::runtime_error("File: " + filename_ + " does not exist!");
    }
  }

  std::ifstream ifp(filename_.c_str());
  char magic_buffer[sizeof(kStackFileMagic)];
  ifp.read(magic_buffer, sizeof(magic_buffer));

  if (strncmp(magic_buffer, kStackFileMagic, sizeof(kStackFileMagic)) != 0) {
    throw std::runtime_error("File: " + filename_ + " is not a StackFile!");
  }

  header_.reset(new internal::StackHeader());
  size_t header_str_size;
  ifp.read(reinterpret_cast<char*>(&header_str_size), sizeof(header_str_size));
  LOGFN_VAR(header_str_size);

  std::vector<char> header_str(header_str_size);
  ifp.read(header_str.data(), header_str_size);

  if (!header_->ParseFromString(
          std::string(header_str.data(), header_str_size))) {
    header_.reset();
    LOG(FATAL) << "Could not parse header from file " << filename_ << std::endl;
  }

  grid_map_.reset(new GridMap(header_->grid_data()));

  size_t num_cells;
  ifp.read(reinterpret_cast<char*>(&num_cells), sizeof(num_cells));
  CHECK_EQ(num_cells, header_->grid_data().num_active_cells());

  for (size_t i = 0; i < num_cells; i++) {
    GridData::Cell* cell = grid_map_->addCell();
    size_t cell_str_size;
    ifp.read(reinterpret_cast<char*>(&cell_str_size), sizeof(cell_str_size));

    std::vector<char> cell_str(cell_str_size);
    ifp.read(cell_str.data(), cell_str_size);
    if (!cell->ParseFromString(std::string(cell_str.data(), cell_str_size))) {
      LOG(FATAL) << "Could not parse grid cell from file " << filename_
                 << std::endl;
    }
  }
  grid_map_->Finalize();

  initialize();
}

StackFile::StackFile(const std::string& filename, const SegyFile& segyfile) {
  SegyOptions opts;
  opts.setTraceHeaderOffsets(segyfile.guessTraceHeaderOffsets());

  createFromSegy(filename, segyfile, opts);

  initialize();
}

StackFile::StackFile(const std::string& filename,
                     const SegyFile& segyfile,
                     const SegyOptions& opts) {
  createFromSegy(filename, segyfile, opts);

  initialize();
}

void StackFile::createFromSegy(const std::string& filename,
                               const SegyFile& segyfile,
                               const SegyOptions& opts) {
  TIMEIT;
  if (!segyfile.is_open() || !(segyfile.open_mode() & std::ios_base::in)) {
    throw std::runtime_error("SegyFile " + segyfile.name() +
                             " not opened for reading!");
  }

  const SegyFile::BinaryHeader& binary_header = segyfile.getBinaryHeader();

  header_.reset(new internal::StackHeader());
  header_->set_version(kStackFileVersion);
  header_->set_description(segyfile.getTextHeader().toString());
  GridData* grid_data = header_->mutable_grid_data();

  grid_data->set_num_samples(
      binary_header.getValueAtOffset<uint16_t>(kSegyOffsetNumSamples));
  grid_data->set_sampling_interval(
      float(
          binary_header.getValueAtOffset<uint16_t>(kSegyOffsetSampleInterval)) /
      1000.0);

  uint16_t unit_code =
      binary_header.getValueAtOffset<uint16_t>(kSegyOffsetMeasurementUnits);
  switch (unit_code) {
    case 1:
      grid_data->set_units(GridData::METERS);
      break;
    case 2:
      grid_data->set_units(GridData::FEET);
      break;
    default:
      grid_data->set_units(GridData::METERS);
  }

  internal::StackHeader::SliceMetadata* inline_metadata =
      header_->mutable_inline_metadata();
  std::string inline_data_file = filename + "_data";
  inline_metadata->set_binary_file(inline_data_file);

  SegyFile::Trace trace;
  uint64_t num_traces_read = 0;
  std::ofstream inline_bin_file(inline_metadata->binary_file().c_str(),
                                std::ios_base::trunc | std::ios_base::binary);

  grid_map_.reset(new GridMap(*grid_data));

  while (segyfile.read(trace)) {
    CHECK_EQ(grid_data->num_samples(), trace.samples().size());
    const SegyFile::Trace::Header& header = trace.header();

    GridData::Cell* grid_cell = grid_map_->addCell();
    grid_cell->set_x_coordinate(
        header.getCoordinateValue(opts.getTraceHeaderOffset(
            SegyFile::Trace::Header::Attribute::X_COORDINATE)));
    grid_cell->set_y_coordinate(
        header.getCoordinateValue(opts.getTraceHeaderOffset(
            SegyFile::Trace::Header::Attribute::Y_COORDINATE)));

    // If 2D or treating as 2D just re-number for a 2D grid.
    if (opts.is2D()) {
      grid_cell->set_inline_number(1);
      grid_cell->set_crossline_number(num_traces_read + 1);
    } else {
      grid_cell->set_inline_number(
          header.getValueAtOffset<int32_t>(opts.getTraceHeaderOffset(
              SegyFile::Trace::Header::Attribute::INLINE_NUMBER)));
      grid_cell->set_crossline_number(
          header.getValueAtOffset<int32_t>(opts.getTraceHeaderOffset(
              SegyFile::Trace::Header::Attribute::CROSSLINE_NUMBER)));
    }

    VLOG(2) << "Cell: " << (*grid_cell) << std::endl;

    inline_bin_file.write(reinterpret_cast<char*>(trace.samples().data()),
                          trace.samples().size() * sizeof(float));

    ++num_traces_read;
    LOG_EVERY_N(INFO, 100000)
        << google::COUNTER << " traces read." << std::endl;
    segyfile.seek(num_traces_read);
  }
  LOG(INFO) << num_traces_read << " traces read" << std::endl;
  inline_bin_file.close();

  grid_map_->RecomputeGrid();

  (*grid_data) = grid_map_->gridData();

  uint64_t num_bins_in_grid = uint64_t(grid_map_->getNumInlines()) *
                              uint64_t(grid_map_->getNumCrosslines());

  if (num_traces_read > num_bins_in_grid) {
    std::ostringstream ostr;
    ostr << "Error: Number of bins in grid is less than the number of traces "
            "read!"
         << std::endl
         << "Could not find valid inline/crossline numbers in file "
         << segyfile.name() << "!" << std::endl
         << "Computed grid:" << std::endl
         << (*grid_data) << std::endl;
    throw std::runtime_error(ostr.str());
  }

  computeInlineMetadata(inline_metadata);

  internal::StackHeader::SliceMetadata* crossline_metadata =
      header_->mutable_crossline_metadata();
  std::string crossline_data_file = filename + "_data_xline";
  crossline_metadata->set_binary_file(crossline_data_file);

  computeCrosslineMetadata(crossline_metadata);

  internal::StackHeader::SliceMetadata* depth_metadata =
      header_->mutable_depth_metadata();
  std::string depth_data_file = filename + "_data_depth";
  depth_metadata->set_binary_file(depth_data_file);

  computeDepthSliceMetadata(depth_metadata);

  internal::UTMZone* utm_zone = header_->mutable_utm_zone();
  utm_zone->set_number(opts.getUtmZone().number());
  std::string letter(1, opts.getUtmZone().letter());
  utm_zone->set_letter(letter.c_str());

  std::ofstream hdr_fp(filename.c_str(),
                       std::ios_base::trunc | std::ios_base::binary);
  hdr_fp.write(kStackFileMagic, sizeof(kStackFileMagic));

  std::string header_str;
  if (!header_->SerializeToString(&header_str)) {
    LOG(FATAL) << "Failed to serialize header for writing" << std::endl;
  }
  size_t header_size = header_str.size();
  hdr_fp.write(reinterpret_cast<char*>(&header_size), sizeof(header_size));
  hdr_fp.write(header_str.c_str(), header_size);

  size_t num_cells = grid_map_->cells().size();
  CHECK_EQ(num_cells, header_->grid_data().num_active_cells());
  hdr_fp.write(reinterpret_cast<char*>(&num_cells), sizeof(num_cells));

  for (const GridData::Cell& cell : grid_map_->cells()) {
    std::string cell_str;
    if (!cell.SerializeToString(&cell_str)) {
      LOG(FATAL) << "Failed to serialize cell for writing" << std::endl;
    }
    size_t cell_str_size = cell_str.size();
    hdr_fp.write(reinterpret_cast<char*>(&cell_str_size),
                 sizeof(cell_str_size));
    hdr_fp.write(cell_str.c_str(), cell_str_size);
  }

  hdr_fp.close();
}

void StackFile::initialize() {
  grid_.reset(new Grid(
      grid_map_.get(),
      UTMZone(header_->utm_zone().number(), header_->utm_zone().letter()[0])));

  data_file_.reset(new MmapFile(header_->inline_metadata().binary_file()));
  if (!data_file_->exists())
    throw std::runtime_error("Missing data file: " + data_file_->name());

  data_file_->open(std::ios_base::in);
  data_file_->map();

  grid_map_->setInlineDataAddress(data_file_->char_addr());

  data_xl_file_.reset(
      new MmapFile(header_->crossline_metadata().binary_file()));
  if (!data_xl_file_->exists())
    data_xl_file_.reset();

  data_ds_file_.reset(new MmapFile(header_->depth_metadata().binary_file()));
  if (!data_ds_file_->exists())
    data_ds_file_.reset();
}

const StackFile::Grid& StackFile::grid() const {
  CHECK_NOTNULL(grid_);
  return (*grid_);
}

void StackFile::readInline(int il,
                           std::vector<float>& data,
                           float fill_value) const {
  readInline(il, data.data(), data.size(), fill_value);
}

void StackFile::readInline(int il,
                           float* buffer,
                           size_t buffer_size,
                           float fill_value) const {
  TIMEIT;
  CHECK_NOTNULL(buffer);

  if (il < grid_->inlineMin() || il > grid_->inlineMax()) {
    throw std::runtime_error("Inline " + std::to_string(il) + " out of range!");
  }

  size_t expected_size = grid_map_->getNumCrosslines() * grid_->numSamples();
  if (buffer_size < expected_size) {
    std::ostringstream ostr;
    ostr << "readInline: data buffer length = " << buffer_size
         << " less than required = " << expected_size;
    throw std::runtime_error(ostr.str());
  }

  std::fill(buffer, buffer + buffer_size, fill_value);

  int offset_index = grid_map_->getActiveInlineIdx(il);
  if (offset_index < 0) {
    return;
  }

  size_t num_trace_bytes = sizeof(float) * grid_->numSamples();
  int num_samples = grid_->numSamples();
  for (int xl = grid_->crosslineMin(), xl_idx = 0; xl <= grid_->crosslineMax();
       xl += grid_->crosslineIncrement(), ++xl_idx) {
    const float* trace = grid_map_->getTrace(il, xl);
    if (trace) {
      ::memcpy(buffer + xl_idx * num_samples, trace, num_trace_bytes);
    }
  }
}

void StackFile::readCrossline(int xl,
                              std::vector<float>& data,
                              float fill_value) const {
  readCrossline(xl, data.data(), data.size(), fill_value);
}

void StackFile::readCrossline(int xl,
                              float* buffer,
                              size_t buffer_size,
                              float fill_value) const {
  TIMEIT;
  CHECK_NOTNULL(buffer);
  if (xl < grid_->crosslineMin() || xl > grid_->crosslineMax()) {
    throw std::runtime_error("Crossline " + std::to_string(xl) +
                             " out of range!");
  }

  size_t expected_size = grid_map_->getNumInlines() * grid_->numSamples();
  if (buffer_size < expected_size) {
    std::ostringstream ostr;
    ostr << "readCrossline: data buffer length = " << buffer_size
         << " less than required = " << expected_size;
    throw std::runtime_error(ostr.str());
  }

  std::fill(buffer, buffer + buffer_size, fill_value);

  if (data_xl_file_) {
    readCrosslineOptimized(xl, buffer, buffer_size);
  } else {
    readCrosslineFromInlineData(xl, buffer, buffer_size);
  }
}

void StackFile::readCrosslineOptimized(int xl,
                                       float* buffer,
                                       size_t buffer_size) const {
  TIMEIT;
  CHECK_NOTNULL(data_xl_file_);
  if (!data_xl_file_->is_open())
    data_xl_file_->open(std::ios_base::in);
  if (!data_xl_file_->is_mapped())
    data_xl_file_->map();

  int offset_index = grid_map_->getActiveCrosslineIdx(xl);
  if (offset_index < 0) {
    return;
  }

  char* trace_addr = data_xl_file_->char_addr() +
                     header_->crossline_metadata().offset(offset_index);
  size_t num_trace_bytes = sizeof(float) * grid_->numSamples();

  for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
       il += grid_->inlineIncrement()) {
    int il_0 = (il - grid_->inlineMin()) / grid_->inlineIncrement();
    if (grid_map_->isCellActive(il, xl)) {
      float* dest_trace = buffer + il_0 * grid_->numSamples();
      CHECK_LE(dest_trace, buffer + buffer_size);
      ::memcpy(dest_trace, trace_addr, num_trace_bytes);
      trace_addr += num_trace_bytes;
    }
  }
}

void StackFile::readCrosslineFromInlineData(int xl,
                                            float* buffer,
                                            size_t buffer_size) const {
  TIMEIT;
  for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
       il += grid_->inlineIncrement()) {
    const float* trace = grid_map_->getTrace(il, xl);
    int il_0 = grid_map_->getInlineIdx(il);
    if (trace) {
      float* dest_trace = buffer + il_0 * grid_->numSamples();
      CHECK_LE(dest_trace, buffer + buffer_size);
      ::memcpy(dest_trace, trace, sizeof(float) * grid_->numSamples());
    }
  }
}

void StackFile::readDepthSlice(unsigned int sample_index,
                               std::vector<float>& data,
                               float fill_value) const {
  readDepthSlice(sample_index, data.data(), data.size(), fill_value);
}

void StackFile::readDepthSlice(unsigned int sample_index,
                               float* buffer,
                               size_t buffer_size,
                               float fill_value) const {
  TIMEIT;
  CHECK_NOTNULL(buffer);
  if (sample_index >= grid_->numSamples()) {
    throw std::runtime_error("Sample index " + std::to_string(sample_index) +
                             " out of range!");
  }

  size_t expected_size =
      grid_map_->getNumInlines() * grid_map_->getNumCrosslines();
  if (buffer_size < expected_size) {
    std::ostringstream ostr;
    ostr << "readDepthSlice: data buffer length = " << buffer_size
         << " less than required = " << expected_size;
    throw std::runtime_error(ostr.str());
  }

  std::fill(buffer, buffer + buffer_size, fill_value);

  if (data_ds_file_) {
    readDepthSliceOptimized(sample_index, buffer, buffer_size);
  } else {
    readDepthSliceFromInlineData(sample_index, buffer, buffer_size);
  }
}

void StackFile::setCrosslineAccessOptimization(bool value) {
  if (value) {
    if (!data_xl_file_) {
      writeCrosslineSlices();
      CHECK_NOTNULL(data_xl_file_);
    }
  } else {
    if (!data_xl_file_)
      return;
    data_xl_file_->remove();
    data_xl_file_.reset();
  }
}

void StackFile::setDepthSliceAccessOptimization(bool value) {
  if (value) {
    if (!data_ds_file_) {
      writeDepthSlices();
      CHECK_NOTNULL(data_ds_file_);
    }
  } else {
    if (!data_ds_file_)
      return;
    data_ds_file_->remove();
    data_ds_file_.reset();
  }
}

void StackFile::writeCrosslineSlices() {
  TIMEIT;
  std::ofstream xl_data_fp(header_->crossline_metadata().binary_file().c_str(),
                           std::ios_base::trunc | std::ios_base::binary);

  size_t num_trace_bytes = sizeof(float) * grid_->numSamples();
  for (int xl = grid_->crosslineMin(); xl <= grid_->crosslineMax();
       xl += grid_->crosslineIncrement()) {
    for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
         il += grid_->inlineIncrement()) {
      const float* trace = grid_map_->getTrace(il, xl);
      if (trace) {
        xl_data_fp.write(reinterpret_cast<const char*>(trace), num_trace_bytes);
      }
    }
  }

  xl_data_fp.close();
  data_xl_file_.reset(
      new MmapFile(header_->crossline_metadata().binary_file()));
}

void StackFile::writeDepthSlices() {
  TIMEIT;
  data_ds_file_.reset(new MmapFile(header_->depth_metadata().binary_file()));
  data_ds_file_->open(std::ios_base::trunc | std::ios_base::binary);

  size_t depth_slice_bytes = 0;
  for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
       il += grid_->inlineIncrement()) {
    for (int xl = grid_->crosslineMin(); xl <= grid_->crosslineMax();
         xl += grid_->crosslineIncrement()) {
      if (grid_map_->isCellActive(il, xl)) {
        depth_slice_bytes += sizeof(float);
      }
    }
  }

  data_ds_file_->expand(depth_slice_bytes * grid_->numSamples());
  data_ds_file_->map();

  std::vector<float> buffer;
  std::vector<float> buffer_ix_trc;
  size_t il_bytes_per_ds_written = 0;
  for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
       il += grid_->inlineIncrement()) {
    int active_il_index = grid_map_->getActiveInlineIdx(il);
    if (active_il_index < 0)
      continue;

    char* il_addr = data_file_->char_addr() +
                    header_->inline_metadata().offset(active_il_index);

    buffer.resize(header_->inline_metadata().size(active_il_index));
    ::memcpy(&buffer[0], il_addr,
             header_->inline_metadata().size(active_il_index));

    size_t trc_num = 0;
    for (size_t iz = 0; iz < grid_->numSamples(); iz++) {
      buffer_ix_trc.resize(grid_map_->getNumCrosslines());
      std::fill(buffer_ix_trc.begin(), buffer_ix_trc.end(), 0.0f);

      trc_num = 0;
      for (int xl = grid_->crosslineMin(); xl <= grid_->crosslineMax();
           xl += grid_->crosslineIncrement()) {
        if (grid_map_->isCellActive(il, xl)) {
          buffer_ix_trc[trc_num] = buffer[trc_num * grid_->numSamples() + iz];
          trc_num++;
        }
      }

      char* ds_addr = data_ds_file_->char_addr() + iz * depth_slice_bytes +
                      il_bytes_per_ds_written;
      ::memcpy(ds_addr, &buffer_ix_trc[0], trc_num * sizeof(float));
    }

    il_bytes_per_ds_written += trc_num * sizeof(float);
  }

  data_ds_file_->unmap();
  data_ds_file_->close();

  data_ds_file_.reset(new MmapFile(header_->depth_metadata().binary_file()));
}

void StackFile::readDepthSliceOptimized(unsigned int sample_index,
                                        float* buffer,
                                        size_t buffer_size) const {
  TIMEIT;
  CHECK_NOTNULL(data_ds_file_);
  if (!data_ds_file_->is_open())
    data_ds_file_->open(std::ios_base::in);
  if (!data_ds_file_->is_mapped())
    data_ds_file_->map();

  std::vector<float> depth_slice(grid_->grid_data_.num_active_cells());
  CHECK_EQ(sizeof(float) * depth_slice.size(),
           size_t(header_->depth_metadata().size(sample_index)));

  ::memcpy(depth_slice.data(),
           data_ds_file_->char_addr() +
               header_->depth_metadata().offset(sample_index),
           sizeof(float) * depth_slice.size());

  size_t src_idx = 0, dest_idx = 0;
  for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
       il += grid_->inlineIncrement()) {
    for (int xl = grid_->crosslineMin(); xl <= grid_->crosslineMax();
         xl += grid_->crosslineIncrement()) {
      if (grid_map_->isCellActive(il, xl)) {
        CHECK_LT(dest_idx, buffer_size);
        buffer[dest_idx] = depth_slice[src_idx++];
      }
      dest_idx++;
    }
  }
}

void StackFile::readDepthSliceFromInlineData(unsigned int sample_index,
                                             float* buffer,
                                             size_t buffer_size) const {
  TIMEIT;
  size_t dest_idx = 0;
  for (int il = grid_->inlineMin(); il <= grid_->inlineMax();
       il += grid_->inlineIncrement()) {
    for (int xl = grid_->crosslineMin(); xl <= grid_->crosslineMax();
         xl += grid_->crosslineIncrement()) {
      if (grid_map_->isCellActive(il, xl)) {
        CHECK_LT(dest_idx, buffer_size);
        const float* trace = grid_map_->getTrace(il, xl);
        buffer[dest_idx] = trace[sample_index];
      }
      dest_idx++;
    }
  }
}

}  // namespace segystack
