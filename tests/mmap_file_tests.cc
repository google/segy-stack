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

#include <fstream>

#include "gtest/gtest.h"
#include "mmap_file.h"
#include "segy_file.h"
#include "test/test_base.h"

namespace segystack {
namespace test {

class MmapFileTest : public ::testing::Test {
  void SetUp() override {
    tmp_file_ = "/tmp/segystack_mmap_file_test.bin";
    std::ofstream ofs(tmp_file_.c_str(),
                      std::ofstream::trunc | std::ofstream::binary);
    test_header_ = "Test data";
    ofs.write(test_header_.c_str(), test_header_.size());
    ofs.close();
  }

  void TearDown() override { ::unlink(tmp_file_.c_str()); }

 protected:
  std::string tmp_file_;
  std::string test_header_;
};

TEST_F(MmapFileTest, NonExistent) {
  std::string nonexistent_file("/tmp/non_existent_file.bin");
  MmapFile fp(nonexistent_file);
  ASSERT_FALSE(fp.exists());
  ASSERT_THROW(fp.open(), std::runtime_error);
  ASSERT_THROW(fp.fd(), std::runtime_error);

  std::ofstream ofp(nonexistent_file.c_str());
  ofp << "data";
  ofp.close();

  ASSERT_NO_THROW(fp.open(std::ios_base::in));
  ASSERT_NO_THROW(fp.close());
  ASSERT_NO_THROW(fp.remove());
  ASSERT_FALSE(fp.exists());
}

TEST_F(MmapFileTest, BasicOperations) {
  MmapFile fp(tmp_file_);
  ASSERT_EQ(fp.name(), tmp_file_);
  ASSERT_TRUE(fp.exists());
  ASSERT_FALSE(fp.is_open());
  ASSERT_FALSE(fp.is_mapped());
  ASSERT_THROW(fp.fd(), std::runtime_error);
  ASSERT_NO_THROW(fp.open(std::ios_base::in));
  ASSERT_TRUE(fp.is_open());
  ASSERT_EQ(fp.size(), test_header_.size());

  struct stat statbuf;
  ASSERT_EQ(::fstat(fp.fd(), &statbuf), 0);

  ASSERT_FALSE(fp.is_mapped());
  ASSERT_NO_THROW(fp.map());
  ASSERT_TRUE(fp.is_mapped());

  ASSERT_NE(fp.addr(), nullptr);
  ASSERT_NE(fp.char_addr(), nullptr);

  std::vector<char> buffer(fp.size());
  ::memcpy(buffer.data(), fp.addr(), fp.size());
  std::string buffer_str(buffer.data(), buffer.size());
  ASSERT_EQ(buffer_str, test_header_);

  ASSERT_NO_THROW(fp.unmap());
  ASSERT_FALSE(fp.is_mapped());
  ASSERT_TRUE(fp.is_open());
  ASSERT_NO_THROW(fp.close());
  ASSERT_FALSE(fp.is_open());
  ASSERT_FALSE(fp.is_mapped());
}

TEST_F(MmapFileTest, CheckTestSegy) {
  int num_samples = 7;
  int num_il = 3, num_xl = 4;
  create_test_segy(tmp_file_, num_samples, 4, num_il, 1, num_xl, 1, 0.0f, 0.0f,
                   1.0f, 1.0f);

  MmapFile fp(tmp_file_);
  ASSERT_NO_THROW(fp.open(std::ios_base::in));
  ASSERT_EQ(
      fp.size(),
      SegyFile::TEXT_HEADER_SIZE + SegyFile::BINARY_HEADER_SIZE +
          num_il * num_xl *
              (SegyFile::TRACE_HEADER_SIZE + sizeof(float) * num_samples));
  ASSERT_NO_THROW(fp.map());
  ASSERT_EQ(fp.end() - fp.begin(), fp.size());
  ASSERT_NO_THROW(fp.close());
  ASSERT_NO_THROW(fp.remove());
  ASSERT_FALSE(fp.exists());
}

}  // namespace test
}  // namespace segystack
