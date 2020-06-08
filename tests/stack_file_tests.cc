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

#include <tuple>

#include "gtest/gtest.h"
#include "stack_file.h"
#include "test/test_base.h"

namespace segystack {
namespace test {

class StackFileTest : public ::testing::Test {
 public:
  void TearDown() override { deleteFile(tmp_file_); }

 protected:
  void runBasicOperationsTest(StackFile& fp) {
    ASSERT_EQ(num_samples_, fp.grid().numSamples());
    ASSERT_FLOAT_EQ(samp_int_ / 1000.0f, fp.grid().samplingInterval());
    ASSERT_EQ(num_il_, fp.grid().numInlines());
    ASSERT_EQ(num_xl_, fp.grid().numCrosslines());

    if (fp.grid().numInlines() > 1) {
      ASSERT_EQ(il_incr_, fp.grid().inlineIncrement());
      ASSERT_FLOAT_EQ(il_spacing_, fp.grid().inlineSpacing());
    }
    if (fp.grid().numCrosslines() > 1) {
      ASSERT_EQ(xl_incr_, fp.grid().crosslineIncrement());
      ASSERT_FLOAT_EQ(xl_spacing_, fp.grid().crosslineSpacing());
    }

    const StackFile::Grid& grid = fp.grid();
    for (int il = grid.inlineMin(); il <= grid.inlineMax();
         il += grid.inlineIncrement()) {
      for (int xl = grid.crosslineMin(); xl <= grid.crosslineMax();
           xl += grid.crosslineIncrement()) {
        ASSERT_TRUE(grid.isBinActive(il, xl));
        StackFile::Grid::Coordinate coord = grid.getCoordinate(il, xl);
        ASSERT_EQ(coord.inline_num, il);
        ASSERT_EQ(coord.crossline_num, xl);
      }
    }

    StackFile::Grid::BoundingBox bbox = grid.boundingBox();

    ASSERT_EQ(bbox.c1.inline_num, grid.inlineMin());
    ASSERT_EQ(bbox.c2.inline_num, grid.inlineMin());
    ASSERT_EQ(bbox.c3.inline_num, grid.inlineMax());
    ASSERT_EQ(bbox.c4.inline_num, grid.inlineMax());

    ASSERT_EQ(bbox.c1.crossline_num, grid.crosslineMin());
    ASSERT_EQ(bbox.c2.crossline_num, grid.crosslineMax());
    ASSERT_EQ(bbox.c3.crossline_num, grid.crosslineMin());
    ASSERT_EQ(bbox.c4.crossline_num, grid.crosslineMax());
  }

  void deleteFile(const std::string& filename) { ::unlink(filename.c_str()); }

  int num_samples_, samp_int_, num_il_, il_incr_, num_xl_, xl_incr_;
  float il_spacing_, xl_spacing_;
  std::string tmp_file_;
};

class StackFileGridTest
    : public StackFileTest,
      public ::testing::WithParamInterface<std::tuple<int, int, int>> {
 public:
  void SetUp() override {
    num_il_ = std::get<0>(GetParam());
    num_xl_ = std::get<1>(GetParam());
    num_samples_ = std::get<2>(GetParam());
    samp_int_ = std::max(num_il_, num_xl_) / std::min(num_il_, num_xl_);
    il_incr_ = std::max(num_il_, num_xl_) + std::min(num_il_, num_xl_);
    xl_incr_ = std::max(num_il_, num_xl_) - std::min(num_il_, num_xl_);
    il_spacing_ = il_incr_ * 20.0;
    xl_spacing_ = xl_incr_ * 10.0;
    std::string suffix = std::to_string(num_il_) + "_" +
                         std::to_string(num_xl_) + "_" +
                         std::to_string(num_samples_);
    tmp_file_ = "/tmp/segystack_stack_file_grid_test" + suffix + ".sgy";
    create_test_segy(tmp_file_, num_samples_, samp_int_, num_il_, il_incr_,
                     num_xl_, xl_incr_, 0.0f, 0.0f, il_spacing_, xl_spacing_);
  }
};

TEST_P(StackFileGridTest, BasicOperations) {
  std::string outfile = tmp_file_ + ".out_stack";
  SegyFile segyfile(tmp_file_);
  segyfile.open(std::ios_base::in);
  StackFile fp(outfile, segyfile, StackFile::SegyOptions());
  segyfile.close();

  runBasicOperationsTest(fp);
  deleteFile(outfile);
  deleteFile(outfile + "_data");
}

INSTANTIATE_TEST_SUITE_P(
    All,
    StackFileGridTest,
    ::testing::Combine(::testing::Values(1, 5, 10, 57, 130),
                       ::testing::Values(1, 4, 11, 61, 84),
                       ::testing::Values(1, 20, 400, 1033)));

class StackFileTraceHeaderTest
    : public StackFileTest,
      public ::testing::WithParamInterface<
          std::tuple<std::pair<int, int>, std::pair<int, int>>> {
 public:
  void SetUp() override {
    int x_offset = std::get<0>(GetParam()).first;
    int y_offset = std::get<0>(GetParam()).second;
    int il_offset = std::get<1>(GetParam()).first;
    int xl_offset = std::get<1>(GetParam()).second;

    StackFile::SegyOptions opts;
    opts.setTraceHeaderOffset(SegyFile::Trace::Header::Attribute::X_COORDINATE,
                              x_offset);
    opts.setTraceHeaderOffset(SegyFile::Trace::Header::Attribute::Y_COORDINATE,
                              y_offset);
    opts.setTraceHeaderOffset(SegyFile::Trace::Header::Attribute::INLINE_NUMBER,
                              il_offset);
    opts.setTraceHeaderOffset(
        SegyFile::Trace::Header::Attribute::CROSSLINE_NUMBER, xl_offset);

    num_il_ = il_offset;
    num_xl_ = xl_offset;
    num_samples_ = x_offset;
    samp_int_ = std::min(num_il_, num_xl_);
    il_incr_ = num_xl_;
    xl_incr_ = num_il_;
    il_spacing_ = il_incr_ * 3.0;
    xl_spacing_ = xl_incr_ * 2.0;
    std::string suffix =
        std::to_string(x_offset) + "_" + std::to_string(y_offset) + "_" +
        std::to_string(il_offset) + "_" + std::to_string(xl_offset);
    tmp_file_ = "/tmp/segystack_stack_file_trc_hdr_test" + suffix + ".sgy";
    create_test_segy(tmp_file_, num_samples_, samp_int_, num_il_, il_incr_,
                     num_xl_, xl_incr_, float(il_incr_), float(xl_incr_),
                     il_spacing_, xl_spacing_, opts);
  }
};

TEST_P(StackFileTraceHeaderTest, BasicOperations) {
  std::string outfile = tmp_file_ + ".out_stack";
  SegyFile segyfile(tmp_file_);
  segyfile.open(std::ios_base::in);
  StackFile fp(outfile, segyfile);
  segyfile.close();

  runBasicOperationsTest(fp);
  deleteFile(outfile);
  deleteFile(outfile + "_data");
}

INSTANTIATE_TEST_SUITE_P(
    All,
    StackFileTraceHeaderTest,
    ::testing::Combine(::testing::Values(std::make_pair(181, 185),
                                         std::make_pair(73, 77),
                                         std::make_pair(201, 205)),
                       ::testing::Values(std::make_pair(189, 193),
                                         std::make_pair(213, 217),
                                         std::make_pair(9, 21))));

}  // namespace test
}  // namespace segystack