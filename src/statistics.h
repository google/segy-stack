/*
 * Copyright 2021 Google LLC
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

#ifndef SEGYSTACK_STATISTICS_H_
#define SEGYSTACK_STATISTICS_H_

#include <string>

namespace segystack {

class Statistics {
 public:
  Statistics(const std::string& stat_file_path)
      : stat_file_path_(stat_file_path) {}
  ~Statistics() {}

  void processInlineTrace(const float* samples, size_t num_elems) {
    processTrace(samples, num_elems, 1);
  }

  void processCrosslineTrace(const float* samples, size_t num_elems) {
    processTrace(samples, num_elems, 2);
  }

  void processDepthTrace(const float* samples, size_t num_elems) {
    processTrace(samples, num_elems, 3);
  }

 private:
  // pass: iline = 1, xline = 2, depth = 3
  void processTrace(const float* samples, size_t num_elems, int pass) {
    for (size_t i = 0; i < num_elems; i++)
      calculateStatistics(samples[i], pass);
  }

  void calculateStatistics(float value, int pass) {}

  std::string stat_file_path_;
};

}  // namespace segystack

#endif  // SEGYSTACK_STATISTICS_H_
