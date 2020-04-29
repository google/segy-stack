#ifndef SEGYSTACK_TESTS_TEST_BASE_H_
#define SEGYSTACK_TESTS_TEST_BASE_H_

#include <string>

#include "stack_file.h"

namespace segystack {
namespace test {

void create_test_segy(
    const std::string& outfile,
    int num_samples,
    int sampling_interval,
    int num_il,
    int il_increment,
    int num_xl,
    int xl_increment,
    const StackFile::SegyOptions& opts = StackFile::SegyOptions());
}
}  // namespace segystack

#endif