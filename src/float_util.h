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

#ifndef SEGYSTACK_FLOAT_UTIL_H_
#define SEGYSTACK_FLOAT_UTIL_H_

#include <cmath>
#include <utility>

namespace segystack {

template <typename T>
T swap_endianness(T value) {
  const int num_bytes = sizeof(T);
  char* bytes = reinterpret_cast<char*>(&value);

  for (int i = 0; i < num_bytes / 2; i++) {
    std::swap(bytes[i], bytes[num_bytes - 1 - i]);
  }
  return value;
}

template <typename T>
T fix_endianness_if_needed(T value) {
#if __BYTE_ORDER == __LITTLE_ENDIAN
  value = swap_endianness(value);
#endif
  return value;
}

// Convert from IBM floating point format to IEEE format.
// IBM: https://en.wikipedia.org/wiki/IBM_hexadecimal_floating_point
// IEEE: https://en.wikipedia.org/wiki/IEEE_754-1985
float ibm_to_ieee(float value, bool is_big_endian_input = true);

}  // namespace segystack

#endif
