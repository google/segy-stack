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

#ifndef SEGYSTACK_LOGGING_H_
#define SEGYSTACK_LOGGING_H_

#include "glog/logging.h"

#include <chrono>

namespace segystack {

namespace internal {
class Timer {
 public:
  Timer(const std::string& name)
      : name_(name), start_time_(std::chrono::system_clock::now()) {
    VLOG(1) << name_ << " called" << std::endl;
  }

  ~Timer() {
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time_;
    VLOG(1) << name_ << " done. Elapsed time: " << elapsed_time.count()
            << " seconds" << std::endl;
  }

 private:
  std::string name_;
  std::chrono::time_point<std::chrono::system_clock> start_time_;
};
}  // namespace internal

}  // namespace segystack

#define TIMEIT segystack::internal::Timer timeit_timer_(__PRETTY_FUNCTION__);

#define LOGFN_VAR(var)                                            \
  {                                                               \
    VLOG(1) << __PRETTY_FUNCTION__ << ":" << #var << " = " << var \
            << std::endl;                                         \
  }

#endif
