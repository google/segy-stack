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

#include "test_base.h"

#include <algorithm>
#include <fstream>

#include "segy_file.h"

namespace segystack {
namespace test {

void create_test_segy(const std::string& outfile,
                      int num_samples,
                      int sampling_interval,
                      int num_il,
                      int il_increment,
                      int num_xl,
                      int xl_increment,
                      float x_origin,
                      float y_origin,
                      float il_spacing,
                      float xl_spacing,
                      const StackFile::SegyOptions& opts) {
  std::ofstream ofs(outfile.c_str(),
                    std::ios_base::trunc | std::ios_base::binary);

  std::string header(3200, ' ');
  std::copy(outfile.begin(), outfile.end(), header.begin());
  ofs.write(header.c_str(), header.size());

  std::vector<char> binary_header(400, 0);

  uint16_t* num_trc_samples =
      reinterpret_cast<uint16_t*>(binary_header.data() + 20);
  (*num_trc_samples) = swap_endianness(static_cast<uint16_t>(num_samples));

  uint16_t* hdr_sampling_interval =
      reinterpret_cast<uint16_t*>(binary_header.data() + 16);
  (*hdr_sampling_interval) =
      swap_endianness(static_cast<uint16_t>(sampling_interval));

  uint16_t* data_format =
      reinterpret_cast<uint16_t*>(binary_header.data() + 24);
  (*data_format) = swap_endianness(uint16_t(5));

  uint8_t* segy_major = reinterpret_cast<uint8_t*>(binary_header.data() + 300);
  (*segy_major) = 1;

  uint16_t* fixed_trc_flag =
      reinterpret_cast<uint16_t*>(binary_header.data() + 302);
  (*fixed_trc_flag) = swap_endianness(uint16_t(1));

  uint16_t* num_ext_hdrs =
      reinterpret_cast<uint16_t*>(binary_header.data() + 304);
  (*num_ext_hdrs) = 0;

  ofs.write(binary_header.data(), binary_header.size());

  std::vector<char> trace_header(240, 0);
  std::vector<float> samples(num_samples, 0);

#define WRITE_VALUE(addr, value)                      \
  {                                                   \
    int32_t* attr = reinterpret_cast<int32_t*>(addr); \
    (*attr) = swap_endianness(int32_t(value));        \
  }

  for (int il = 0; il < num_il; ++il) {
    for (int xl = 0; xl < num_xl; ++xl) {
      char* il_num_addr =
          trace_header.data() +
          opts.getTraceHeaderOffset(
              SegyFile::Trace::Header::Attribute::INLINE_NUMBER) -
          1;
      WRITE_VALUE(il_num_addr, (il * il_increment));

      char* xl_num_addr =
          trace_header.data() +
          opts.getTraceHeaderOffset(
              SegyFile::Trace::Header::Attribute::CROSSLINE_NUMBER) -
          1;
      WRITE_VALUE(xl_num_addr, (xl * xl_increment));

      char* x_coord_addr =
          trace_header.data() +
          opts.getTraceHeaderOffset(
              SegyFile::Trace::Header::Attribute::X_COORDINATE) -
          1;
      WRITE_VALUE(x_coord_addr, (x_origin + xl * xl_spacing));

      char* y_coord_addr =
          trace_header.data() +
          opts.getTraceHeaderOffset(
              SegyFile::Trace::Header::Attribute::Y_COORDINATE) -
          1;
      WRITE_VALUE(y_coord_addr, (y_origin + il * il_spacing));

      ofs.write(trace_header.data(), trace_header.size());
      ofs.write(reinterpret_cast<char*>(samples.data()),
                sizeof(float) * samples.size());
    }
  }

  ofs.close();
}

}  // namespace test
}  // namespace segystack
