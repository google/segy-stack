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

#include "float_util.h"

#include <cstdint>
#include <limits>

#define SIGN_BIT_MASK 0x80000000
#define EXPONENT_BIT_MASK 0x7f000000
#define FRACTION_BIT_MASK 0x00ffffff
#define FRACTION_IMPLICIT_BIT_MASK 0x00800000

namespace segystack {

float ibm_to_ieee(float value, bool is_big_endian_input) {
  if (is_big_endian_input) {
    value = swap_endianness(value);
  }

  int32_t* int_addr = reinterpret_cast<int32_t*>(&value);
  int32_t int_val = *int_addr;

  int32_t sign = int_val >> 31;
  int32_t fraction = int_val & FRACTION_BIT_MASK;

  if (fraction == 0) {
    return sign ? -0.0f : 0.0f;
  }

  // Convert exponent to be of base 2 and remove IBM exponent bias.
  int32_t exponent = ((int_val & EXPONENT_BIT_MASK) >> 22) - 256;

  // Drop the last bit since we can store only 23 bits in IEEE.
  fraction >>= 1;

  // Normalize such that the implicit leading bit of the fraction is 1.
  while (fraction && (fraction & FRACTION_IMPLICIT_BIT_MASK) == 0) {
    fraction <<= 1;
    --exponent;
  }

  // Drop the implicit leading bit.
  fraction &= 0x007fffff;

  // Add IEEE bias to the exponent.
  exponent += 127;

  // Handle overflow.
  if (exponent >= 255) {
    return (sign ? -std::numeric_limits<float>::max()
                 : std::numeric_limits<float>::max());
  }

  int32_t ieee_value;

  // Handle underflow.
  if (exponent <= 0)
    ieee_value = (sign << 31) | fraction;
  else
    ieee_value = (sign << 31) | (exponent << 23) | fraction;

  float* float_addr = reinterpret_cast<float*>(&ieee_value);
  return *float_addr;
}

}  // namespace segystack
