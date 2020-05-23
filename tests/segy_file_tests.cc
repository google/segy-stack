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
#include "segy_file.h"
#include "test/test_base.h"

namespace segystack {
namespace test {

class SegyFileTest
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<int, int, int>> {
 public:
  void SetUp() override {
    num_il_ = std::get<0>(GetParam());
    num_xl_ = std::get<1>(GetParam());
    num_samples_ = std::get<2>(GetParam());
    samp_int_ = std::max(num_il_, num_xl_);
    std::string suffix = std::to_string(num_il_) + "_" +
                         std::to_string(num_xl_) + "_" +
                         std::to_string(num_samples_);
    tmp_file_ = "/tmp/segystack_segy_file_test" + suffix + ".sgy";
    create_test_segy(tmp_file_, num_samples_, samp_int_, num_il_, 1, num_xl_, 1,
                     0.0f, 0.0f, 1.0f, 1.0f);
  }

  void TearDown() override { ::unlink(tmp_file_.c_str()); }

 protected:
  int num_samples_, samp_int_, num_il_, num_xl_;
  std::string tmp_file_;
};

TEST_P(SegyFileTest, SanityChecks) {
  ASSERT_EQ(sizeof(SegyFile::TextHeader), SegyFile::TEXT_HEADER_SIZE);
  ASSERT_EQ(sizeof(SegyFile::BinaryHeader), SegyFile::BINARY_HEADER_SIZE);
  ASSERT_EQ(sizeof(SegyFile::Trace::Header), SegyFile::TRACE_HEADER_SIZE);
}

TEST_P(SegyFileTest, BasicOperations) {
  SegyFile fp(tmp_file_);
  ASSERT_NO_THROW(fp.open(std::ios_base::in));
  ASSERT_NE(fp.getTextHeader().toString().find(tmp_file_), std::string::npos);
  ASSERT_EQ(fp.getBinaryHeader().getValueAtOffset<uint16_t>(21), num_samples_);
  ASSERT_EQ(fp.getBinaryHeader().getValueAtOffset<uint16_t>(17), samp_int_);

  SegyFile::Trace trace;
  int num_trcs_read = 0;
  while (fp.read(trace)) {
    ++num_trcs_read;
    ASSERT_EQ(trace.samples().size(), num_samples_);
    ASSERT_NO_THROW(fp.seek(num_trcs_read));
  }
  ASSERT_EQ(num_trcs_read, num_il_ * num_xl_);
  ASSERT_NO_THROW(fp.close());
}

INSTANTIATE_TEST_SUITE_P(
    All,
    SegyFileTest,
    ::testing::Combine(::testing::Values(1, 5, 10, 57, 130),
                       ::testing::Values(1, 4, 11, 61, 84),
                       ::testing::Values(1, 20, 400, 1033)));

}  // namespace test
}  // namespace segystack