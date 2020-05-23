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

#include <cmath>
#include <functional>
#include <set>
#include <unordered_map>

#include "logging.h"
#include "segy_file.h"

std::ostream& operator<<(std::ostream& os, const segystack::Grid& grid) {
  os << "IL range : " << grid.inline_min() << " - " << grid.inline_max()
     << std::endl;
  os << "XL range : " << grid.crossline_min() << " - " << grid.crossline_max()
     << std::endl;
  os << "IL increment : " << grid.inline_increment() << std::endl;
  os << "XL increment : " << grid.crossline_increment() << std::endl;
  os << "IL spacing : " << grid.inline_spacing() << std::endl;
  os << "XL spacing : " << grid.crossline_spacing() << std::endl;
  os << "Sampling interval : " << grid.sampling_interval() << std::endl;
  os << "Num depth samples : " << grid.num_samples() << std::endl;
  os << "Units : ";
  switch (grid.units()) {
    case segystack::Grid::METERS:
      os << "Meters";
      break;
    case segystack::Grid::FEET:
      os << "Feet";
      break;
    default:
      os << "Unknown";
  };
  os << std::endl;
  os << "Num active cells : " << grid.num_active_cells() << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const segystack::Grid::Cell& cell) {
  os << "Coordinate (x, y) : (" << cell.x_coordinate() << ", "
     << cell.y_coordinate() << ")" << std::endl;
  os << "Grid numbers (IL, XL) : (" << cell.inline_number() << ", "
     << cell.crossline_number() << ")" << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const segystack::UTMZone& utm) {
  os << utm.number() << " " << utm.letter();
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const segystack::StackFile::SegyOptions& opts) {
  using SegyOptions = segystack::StackFile::SegyOptions;
  os << "UTM zone: " << opts.getUtmZone() << std::endl;
  os << "Inline number offset: "
     << opts.getTraceHeaderOffset(SegyOptions::INLINE_NUMBER) << std::endl;
  os << "Crossline number offset: "
     << opts.getTraceHeaderOffset(SegyOptions::CROSSLINE_NUMBER) << std::endl;
  os << "X Coordinate offset: "
     << opts.getTraceHeaderOffset(SegyOptions::X_COORDINATE) << std::endl;
  os << "Y Coordinate offset: "
     << opts.getTraceHeaderOffset(SegyOptions::Y_COORDINATE) << std::endl;
  return os;
}

namespace segystack {

constexpr char kStackFileMagic[] = "**StackFile**\n";
constexpr int kStackFileVersion = 1;
constexpr int kSegyOffsetSampleInterval = 17;
constexpr int kSegyOffsetNumSamples = 21;
constexpr int kSegyOffsetMeasurementUnits = 55;

StackFile::~StackFile() = default;

StackFile::SegyOptions::SegyOptions() {
  utm_zone_.set_number(32);
  utm_zone_.set_letter("U");

  offsets_[X_COORDINATE] = 181;
  offsets_[Y_COORDINATE] = 185;
  offsets_[INLINE_NUMBER] = 189;
  offsets_[CROSSLINE_NUMBER] = 193;
}

void StackFile::SegyOptions::setUtmZone(int num, char zone) {
  if (num < 1 || num > 60)
    throw std::runtime_error("UTM zone number should be in range [1, 60]");

  std::string letter(1, std::toupper(zone));
  std::set<char> zones;
  const std::set<char> forbidden = {'A', 'B', 'Y', 'Z', 'I', 'O'};
  for (char c = 'A'; c <= 'Z'; c++) {
    if (forbidden.find(c) == forbidden.end())
      zones.insert(c);
  }

  if (zones.find(letter[0]) == zones.end()) {
    std::ostringstream ostr;
    for (char c : zones)
      ostr << c << ", ";
    throw std::runtime_error("UTM zone letter should be one of " + ostr.str());
  }

  utm_zone_.set_number(num);
  utm_zone_.set_letter(letter.c_str());
}

void StackFile::SegyOptions::setTraceHeaderOffset(TraceHeaderAttribute attr,
                                                  int offset) {
  if (offset < 1 || offset > 231)
    throw std::runtime_error("Offset specified outside of range [1, 231]");
  offsets_[attr] = offset;
}

class StackFile::GridMap {
 public:
  GridMap() {}

  GridMap(const Grid& grid) : grid_(grid) {}

  Grid::Cell* addCell() {
    cells_.push_back(Grid::Cell());
    return &(cells_.back());
  }

  void RecomputeGrid() {
    grid_.set_inline_min(std::numeric_limits<int>::max());
    grid_.set_crossline_min(std::numeric_limits<int>::max());
    grid_.set_inline_max(std::numeric_limits<int>::min());
    grid_.set_crossline_max(std::numeric_limits<int>::min());

    for (const Grid::Cell& cell : cells_) {
      grid_.set_inline_min(std::min(grid_.inline_min(), cell.inline_number()));
      grid_.set_crossline_min(
          std::min(grid_.crossline_min(), cell.crossline_number()));
      grid_.set_inline_max(std::max(grid_.inline_max(), cell.inline_number()));
      grid_.set_crossline_max(
          std::max(grid_.crossline_max(), cell.crossline_number()));
    }
    grid_.set_num_active_cells(cells_.size());

    grid_.set_inline_increment(computeIncrement(
        grid_.inline_min(), grid_.inline_max(),
        [](const Grid::Cell& cell) { return cell.inline_number(); }));

    grid_.set_crossline_increment(computeIncrement(
        grid_.crossline_min(), grid_.crossline_max(),
        [](const Grid::Cell& cell) { return cell.crossline_number(); }));

    Finalize();

    grid_.set_inline_spacing(computeSpacing(&GridMap::getCoordInNextInline));
    grid_.set_crossline_spacing(
        computeSpacing(&GridMap::getCoordInNextCrossline));
  }

  void Finalize() {
    cell_map_.resize(getNumInlines());
    for (int i = 0; i < getNumInlines(); i++) {
      cell_map_[i].resize(getNumCrosslines(), nullptr);
    }

    for (const Grid::Cell& cell : cells_) {
      int il_idx = getInlineIdx(cell.inline_number());
      int xl_idx = getCrosslineIdx(cell.crossline_number());
      cell_map_[il_idx][xl_idx] = &cell;
    }

    active_il_index_map_.clear();
    int il_idx = 0;
    for (int il = grid_.inline_min(); il <= grid_.inline_max();
         il += grid_.inline_increment()) {
      for (int xl = grid_.crossline_min(); xl <= grid_.crossline_max();
           xl += grid_.crossline_increment()) {
        if (isCellActive(il, xl)) {
          active_il_index_map_[il] = il_idx++;
          break;
        }
      }
    }

    active_xl_index_map_.clear();
    int xl_idx = 0;
    for (int xl = grid_.crossline_min(); xl <= grid_.crossline_max();
         xl += grid_.crossline_increment()) {
      for (int il = grid_.inline_min(); il <= grid_.inline_max();
           il += grid_.inline_increment()) {
        if (isCellActive(il, xl)) {
          active_xl_index_map_[xl] = xl_idx++;
          break;
        }
      }
    }
  }

  int getNumInlines() const {
    CHECK_GT(grid_.inline_increment(), 0);
    return ((grid_.inline_max() - grid_.inline_min()) /
            grid_.inline_increment()) +
           1;
  }

  int getNumCrosslines() const {
    CHECK_GT(grid_.crossline_increment(), 0);
    return ((grid_.crossline_max() - grid_.crossline_min()) /
            grid_.crossline_increment()) +
           1;
  }

  int getInlineIncrement() const { return grid_.inline_increment(); }
  int getCrosslineIncrement() const { return grid_.crossline_increment(); }

  const Grid& grid() const { return grid_; }

  const std::vector<Grid::Cell>& cells() const { return cells_; }

  bool isCellActive(int il, int xl) const {
    try {
      return cell_map_.at(getInlineIdx(il)).at(getCrosslineIdx(xl)) != nullptr;
    } catch (const std::out_of_range& e) {
      DLOG(INFO) << __FUNCTION__ << ": (" << il << ", " << xl
                 << ") out of bounds";
    }
    return false;
  }

  void setInlineDataAddress(const char* trc_addr) {
    size_t num_trace_bytes = sizeof(float) * grid_.num_samples();

    trace_map_.resize(getNumInlines());
    for (size_t i = 0; i < trace_map_.size(); i++) {
      trace_map_[i].resize(getNumCrosslines(), nullptr);
    }

    for (int il = grid_.inline_min(); il <= grid_.inline_max();
         il += grid_.inline_increment()) {
      for (int xl = grid_.crossline_min(); xl <= grid_.crossline_max();
           xl += grid_.crossline_increment()) {
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
    return (il - grid_.inline_min()) / grid_.inline_increment();
  }

  int getActiveInlineIdx(int il) const {
    if (active_il_index_map_.find(il) != active_il_index_map_.end())
      return active_il_index_map_.at(il);
    return -1;
  }

  int getCrosslineIdx(int xl) const {
    return (xl - grid_.crossline_min()) / grid_.crossline_increment();
  }

  int getActiveCrosslineIdx(int xl) const {
    if (active_xl_index_map_.find(xl) != active_xl_index_map_.end())
      return active_xl_index_map_.at(xl);
    return -1;
  }

 private:
  int computeIncrement(int min_value,
                       int max_value,
                       std::function<int(const Grid::Cell&)> line_number) {
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

  typedef const Grid::Cell*(GridMap::*GetNextCoordMethod)(size_t, size_t) const;
  float computeSpacing(GetNextCoordMethod get_next_coord) {
    float min_spacing = std::numeric_limits<float>::max();
    for (size_t i = 0; i < cell_map_.size(); i++) {
      for (size_t j = 0; j < cell_map_[i].size(); j++) {
        const Grid::Cell* c1 = cell_map_[i][j];
        const Grid::Cell* c2 = (this->*get_next_coord)(i, j);
        if (c1 && c2) {
          float dist =
              std::sqrt(std::pow(c1->x_coordinate() - c2->x_coordinate(), 2.0) +
                        std::pow(c1->y_coordinate() - c2->y_coordinate(), 2.0));
          min_spacing = std::min(min_spacing, dist);
        }
      }
    }

    if (min_spacing >= std::numeric_limits<float>::max()) {
      return 1.0;
    }
    return min_spacing;
  }

  const Grid::Cell* getCoordInNextInline(size_t il_idx, size_t xl_idx) const {
    if (il_idx >= 0 && il_idx < cell_map_.size() - 1 && xl_idx >= 0 &&
        xl_idx < cell_map_[il_idx + 1].size()) {
      return cell_map_[il_idx + 1][xl_idx];
    }
    return nullptr;
  }

  const Grid::Cell* getCoordInNextCrossline(size_t il_idx,
                                            size_t xl_idx) const {
    if (il_idx >= 0 && il_idx < cell_map_.size() && xl_idx >= 0 &&
        xl_idx < cell_map_[il_idx].size() - 1) {
      return cell_map_[il_idx][xl_idx + 1];
    }
    return nullptr;
  }

  Grid grid_;
  std::vector<Grid::Cell> cells_;
  std::vector<std::vector<const Grid::Cell*>> cell_map_;
  std::unordered_map<int, int> active_il_index_map_, active_xl_index_map_;
  std::vector<std::vector<const float*>> trace_map_;
};

void StackFile::computeInlineMetadata(StackHeader::SliceMetadata* il_metadata) {
  CHECK_NOTNULL(grid_map_);

  size_t num_bytes_written = 0;
  const Grid& grid = grid_map_->grid();

  for (int il = grid.inline_min(); il <= grid.inline_max();
       il += grid.inline_increment()) {
    size_t il_size = 0;

    for (int xl = grid.crossline_min(); xl <= grid.crossline_max();
         xl += grid.crossline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        il_size += sizeof(float) * grid.num_samples();
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
    StackHeader::SliceMetadata* xl_metadata) {
  CHECK_NOTNULL(grid_map_);

  size_t num_bytes_written = 0;
  const Grid& grid = grid_map_->grid();

  for (int xl = grid.crossline_min(); xl <= grid.crossline_max();
       xl += grid.crossline_increment()) {
    size_t xl_size = 0;

    for (int il = grid.inline_min(); il <= grid.inline_max();
         il += grid.inline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        xl_size += sizeof(float) * grid.num_samples();
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
    StackHeader::SliceMetadata* depth_metadata) {
  CHECK_NOTNULL(grid_map_);

  size_t depth_slice_bytes = 0;
  const Grid& grid = grid_map_->grid();

  for (int il = grid.inline_min(); il <= grid.inline_max();
       il += grid.inline_increment()) {
    for (int xl = grid.crossline_min(); xl <= grid.crossline_max();
         xl += grid.crossline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        depth_slice_bytes += sizeof(float);
      }
    }
  }

  size_t num_bytes_written = 0;
  for (size_t iz = 0; iz < grid.num_samples(); iz++) {
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

  header_.reset(new StackHeader());
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

  grid_map_.reset(new GridMap(header_->grid()));

  size_t num_cells;
  ifp.read(reinterpret_cast<char*>(&num_cells), sizeof(num_cells));
  CHECK_EQ(num_cells, header_->grid().num_active_cells());

  for (size_t i = 0; i < num_cells; i++) {
    Grid::Cell* cell = grid_map_->addCell();
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

StackFile::StackFile(const std::string& filename,
                     const SegyFile& segyfile,
                     const SegyOptions& opts) {
  TIMEIT;
  if (!segyfile.is_open() || !(segyfile.open_mode() & std::ios_base::in)) {
    throw std::runtime_error("SegyFile " + segyfile.name() + " not opened for reading!");
  }

  const SegyFile::BinaryHeader& binary_header = segyfile.getBinaryHeader();

  header_.reset(new StackHeader());
  header_->set_version(kStackFileVersion);
  header_->set_description(segyfile.getTextHeader().toString());
  Grid* grid = header_->mutable_grid();

  grid->set_num_samples(
      binary_header.getValueAtOffset<uint16_t>(kSegyOffsetNumSamples));
  grid->set_sampling_interval(float(binary_header.getValueAtOffset<uint16_t>(
                                  kSegyOffsetSampleInterval)) /
                              1000.0);

  uint16_t unit_code =
      binary_header.getValueAtOffset<uint16_t>(kSegyOffsetMeasurementUnits);
  switch (unit_code) {
    case 1:
      grid->set_units(Grid_Units_METERS);
      break;
    case 2:
      grid->set_units(Grid_Units_FEET);
      break;
    default:
      grid->set_units(Grid_Units_METERS);
  }

  StackHeader::SliceMetadata* inline_metadata =
      header_->mutable_inline_metadata();
  std::string inline_data_file = filename + "_data";
  inline_metadata->set_binary_file(inline_data_file);

  SegyFile::Trace trace;
  uint64_t num_traces_read = 0;
  std::ofstream inline_bin_file(inline_metadata->binary_file().c_str(),
                                std::ios_base::trunc | std::ios_base::binary);

  grid_map_.reset(new GridMap(*grid));

  while (segyfile.read(trace)) {
    CHECK_EQ(grid->num_samples(), trace.samples().size());
    const SegyFile::Trace::Header& header = trace.header();

    Grid::Cell* grid_cell = grid_map_->addCell();
    grid_cell->set_x_coordinate(header.getCoordinateValue(
        opts.getTraceHeaderOffset(SegyOptions::X_COORDINATE)));
    grid_cell->set_y_coordinate(header.getCoordinateValue(
        opts.getTraceHeaderOffset(SegyOptions::Y_COORDINATE)));
    grid_cell->set_inline_number(header.getValueAtOffset<int32_t>(
        opts.getTraceHeaderOffset(SegyOptions::INLINE_NUMBER)));
    grid_cell->set_crossline_number(header.getValueAtOffset<int32_t>(
        opts.getTraceHeaderOffset(SegyOptions::CROSSLINE_NUMBER)));

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
  (*grid) = grid_map_->grid();

  computeInlineMetadata(inline_metadata);

  StackHeader::SliceMetadata* crossline_metadata =
      header_->mutable_crossline_metadata();
  std::string crossline_data_file = filename + "_data_xline";
  crossline_metadata->set_binary_file(crossline_data_file);

  computeCrosslineMetadata(crossline_metadata);

  StackHeader::SliceMetadata* depth_metadata =
      header_->mutable_depth_metadata();
  std::string depth_data_file = filename + "_data_depth";
  depth_metadata->set_binary_file(depth_data_file);

  computeDepthSliceMetadata(depth_metadata);

  UTMZone* utm_zone = header_->mutable_utm_zone();
  (*utm_zone) = opts.getUtmZone();

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
  CHECK_EQ(num_cells, header_->grid().num_active_cells());
  hdr_fp.write(reinterpret_cast<char*>(&num_cells), sizeof(num_cells));

  for (const Grid::Cell& cell : grid_map_->cells()) {
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

  initialize();
}

void StackFile::initialize() {
  grid_ = &(grid_map_->grid());

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

const Grid& StackFile::grid() const {
  CHECK_NOTNULL(grid_);
  return *grid_;
}

int StackFile::getNumInlines() const {
  return grid_map_->getNumInlines();
}

int StackFile::getNumCrosslines() const {
  return grid_map_->getNumCrosslines();
}

const UTMZone& StackFile::getUtmZone() const {
  CHECK_NOTNULL(header_);
  return header_->utm_zone();
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

  if (il < grid_->inline_min() || il > grid_->inline_max()) {
    throw std::runtime_error("Inline " + std::to_string(il) + " out of range!");
  }

  size_t expected_size = grid_map_->getNumCrosslines() * grid_->num_samples();
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

  size_t num_trace_bytes = sizeof(float) * grid_map_->grid().num_samples();
  int num_samples = grid_map_->grid().num_samples();
  for (int xl = grid_->crossline_min(), xl_idx = 0;
       xl <= grid_->crossline_max();
       xl += grid_->crossline_increment(), ++xl_idx) {
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
  if (xl < grid_->crossline_min() || xl > grid_->crossline_max()) {
    throw std::runtime_error("Crossline " + std::to_string(xl) +
                             " out of range!");
  }

  size_t expected_size = grid_map_->getNumInlines() * grid_->num_samples();
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
  size_t num_trace_bytes = sizeof(float) * grid_map_->grid().num_samples();

  for (int il = grid_->inline_min(); il <= grid_->inline_max();
       il += grid_->inline_increment()) {
    int il_0 = (il - grid_->inline_min()) / grid_->inline_increment();
    if (grid_map_->isCellActive(il, xl)) {
      float* dest_trace = buffer + il_0 * grid_->num_samples();
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
  for (int il = grid_->inline_min(); il <= grid_->inline_max();
       il += grid_->inline_increment()) {
    const float* trace = grid_map_->getTrace(il, xl);
    int il_0 = grid_map_->getInlineIdx(il);
    if (trace) {
      float* dest_trace = buffer + il_0 * grid_->num_samples();
      CHECK_LE(dest_trace, buffer + buffer_size);
      ::memcpy(dest_trace, trace, sizeof(float) * grid_->num_samples());
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
  if (sample_index >= grid_->num_samples()) {
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

  size_t num_trace_bytes = sizeof(float) * grid_->num_samples();
  for (int xl = grid_->crossline_min(); xl <= grid_->crossline_max();
       xl += grid_->crossline_increment()) {
    for (int il = grid_->inline_min(); il <= grid_->inline_max();
         il += grid_->inline_increment()) {
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
  for (int il = grid_->inline_min(); il <= grid_->inline_max();
       il += grid_->inline_increment()) {
    for (int xl = grid_->crossline_min(); xl <= grid_->crossline_max();
         xl += grid_->crossline_increment()) {
      if (grid_map_->isCellActive(il, xl)) {
        depth_slice_bytes += sizeof(float);
      }
    }
  }

  data_ds_file_->expand(depth_slice_bytes * grid_->num_samples());
  data_ds_file_->map();

  std::vector<float> buffer;
  std::vector<float> buffer_ix_trc;
  size_t il_bytes_per_ds_written = 0;
  for (int il = grid_->inline_min(); il <= grid_->inline_max();
       il += grid_->inline_increment()) {
    int active_il_index = grid_map_->getActiveInlineIdx(il);
    if (active_il_index < 0)
      continue;

    char* il_addr = data_file_->char_addr() +
                    header_->inline_metadata().offset(active_il_index);

    buffer.resize(header_->inline_metadata().size(active_il_index));
    ::memcpy(&buffer[0], il_addr,
                header_->inline_metadata().size(active_il_index));

    size_t trc_num = 0;
    for (size_t iz = 0; iz < grid_->num_samples(); iz++) {
      buffer_ix_trc.resize(grid_map_->getNumCrosslines());
      std::fill(buffer_ix_trc.begin(), buffer_ix_trc.end(), 0.0f);

      trc_num = 0;
      for (int xl = grid_->crossline_min(); xl <= grid_->crossline_max();
           xl += grid_->crossline_increment()) {
        if (grid_map_->isCellActive(il, xl)) {
          buffer_ix_trc[trc_num] = buffer[trc_num * grid_->num_samples() + iz];
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

  std::vector<float> depth_slice(grid_->num_active_cells());
  CHECK_EQ(sizeof(float) * depth_slice.size(),
           size_t(header_->depth_metadata().size(sample_index)));

  ::memcpy(depth_slice.data(),
           data_ds_file_->char_addr() +
               header_->depth_metadata().offset(sample_index),
           sizeof(float) * depth_slice.size());

  size_t src_idx = 0, dest_idx = 0;
  for (int il = grid_->inline_min(); il <= grid_->inline_max();
       il += grid_->inline_increment()) {
    for (int xl = grid_->crossline_min(); xl <= grid_->crossline_max();
         xl += grid_->crossline_increment()) {
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
  for (int il = grid_->inline_min(); il <= grid_->inline_max();
       il += grid_->inline_increment()) {
    for (int xl = grid_->crossline_min(); xl <= grid_->crossline_max();
         xl += grid_->crossline_increment()) {
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
