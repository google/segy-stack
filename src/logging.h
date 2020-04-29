#ifndef SEGYSTACK_LOGGING_H_
#define SEGYSTACK_LOGGING_H_

#include <glog/logging.h>

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