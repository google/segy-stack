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
#include <limits>

#include "float_util.h"
#include "gtest/gtest.h"
#include "test/test_base.h"

namespace segystack {
namespace test {

class FloatUtilTest : public ::testing::Test {
 protected:
  float getFloat(int32_t int_val) {
    float* input = reinterpret_cast<float*>(&int_val);
    return *input;
  }
};

TEST_F(FloatUtilTest, Basic) {
  ASSERT_FLOAT_EQ(ibm_to_ieee(getFloat(0x14ad6c42)), 108.6761);
  ASSERT_FLOAT_EQ(ibm_to_ieee(getFloat(0x00a076c2)), -118.625);
  ASSERT_FLOAT_EQ(ibm_to_ieee(getFloat(0xc276a000), false), -118.625);
}

}  // namespace test
}  // namespace segystack
