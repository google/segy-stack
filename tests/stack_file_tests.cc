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

class StackFileTest
    : public ::testing::Test,
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
    tmp_file_ = "/tmp/segystack_stack_file_test" + suffix + ".sgy";
    create_test_segy(tmp_file_, num_samples_, samp_int_, num_il_, il_incr_,
                     num_xl_, xl_incr_, 0.0f, 0.0f, il_spacing_, xl_spacing_);
  }

  void TearDown() override { deleteFile(tmp_file_); }

 protected:
  void deleteFile(const std::string& filename) { ::unlink(filename.c_str()); }

  int num_samples_, samp_int_, num_il_, il_incr_, num_xl_, xl_incr_;
  float il_spacing_, xl_spacing_;
  std::string tmp_file_;
};

TEST_P(StackFileTest, BasicOperations) {
  std::string outfile = tmp_file_ + ".out_stack";
  SegyFile segyfile(tmp_file_);
  segyfile.open(std::ios_base::in);
  StackFile fp(outfile, segyfile, StackFile::SegyOptions());
  segyfile.close();

  ASSERT_EQ(num_samples_, fp.grid().num_samples());
  ASSERT_FLOAT_EQ(samp_int_ / 1000.0f, fp.grid().sampling_interval());
  ASSERT_EQ(num_il_, fp.getNumInlines());
  ASSERT_EQ(num_xl_, fp.getNumCrosslines());

  if (fp.getNumInlines() > 1) {
    ASSERT_EQ(il_incr_, fp.grid().inline_increment());
    ASSERT_FLOAT_EQ(il_spacing_, fp.grid().inline_spacing());
  }
  if (fp.getNumCrosslines() > 1) {
    ASSERT_EQ(xl_incr_, fp.grid().crossline_increment());
    ASSERT_FLOAT_EQ(xl_spacing_, fp.grid().crossline_spacing());
  }

  deleteFile(outfile);
  deleteFile(outfile + "_data");
}

INSTANTIATE_TEST_SUITE_P(
    All,
    StackFileTest,
    ::testing::Combine(::testing::Values(1, 5, 10, 57, 130),
                       ::testing::Values(1, 4, 11, 61, 84),
                       ::testing::Values(1, 20, 400, 1033)));

}  // namespace test
}  // namespace segystack