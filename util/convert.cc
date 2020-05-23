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


#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/flags/usage.h>
#include <glog/logging.h>

#include <iostream>

#include "segy_file.h"
#include "stack_file.h"

ABSL_FLAG(std::string, input_file, "", "Input SEGY file");
ABSL_FLAG(std::string, output_file, "", "Output Stack file");
ABSL_FLAG(int, il_offset, 189, "Inline number byte offset");
ABSL_FLAG(int, xl_offset, 193, "Crossline number byte offset");
ABSL_FLAG(int, x_coord_offset, 181, "X coordinate byte offset");
ABSL_FLAG(int, y_coord_offset, 185, "Y coordinate byte offset");
ABSL_FLAG(int, utm_zone_num, 32, "UTM Zone number");
ABSL_FLAG(std::string, utm_zone_letter, "U", "UTM Zone letter");
ABSL_FLAG(bool, enable_crossline_opt, false, "Enable fast crossline access");
ABSL_FLAG(bool, enable_depth_opt, false, "Enable fast depth slice access");

using namespace segystack;

void usage() {
  std::cout << absl::ProgramUsageMessage() << std::endl;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();

  absl::SetProgramUsageMessage(
      absl::StrCat("Usage:\n", argv[0],
                   " --input_file=segy_file.sgy --output_file=outfile.stack"));
  absl::ParseCommandLine(argc, argv);

  std::string infile;
  if (FLAGS_input_file.IsSpecifiedOnCommandLine()) {
    infile = FLAGS_input_file.Get();
  } else {
    usage();
    return -1;
  }

  std::string outfile = "out.stack";
  if (FLAGS_output_file.IsSpecifiedOnCommandLine())
    outfile = FLAGS_output_file.Get();

  StackFile::SegyOptions opts;
  if (FLAGS_il_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(StackFile::SegyOptions::INLINE_NUMBER,
                              FLAGS_il_offset.Get());
  }

  if (FLAGS_xl_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(StackFile::SegyOptions::CROSSLINE_NUMBER,
                              FLAGS_xl_offset.Get());
  }

  if (FLAGS_x_coord_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(StackFile::SegyOptions::X_COORDINATE,
                              FLAGS_x_coord_offset.Get());
  }

  if (FLAGS_y_coord_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(StackFile::SegyOptions::Y_COORDINATE,
                              FLAGS_y_coord_offset.Get());
  }

  if (FLAGS_utm_zone_num.IsSpecifiedOnCommandLine() &&
      FLAGS_utm_zone_letter.IsSpecifiedOnCommandLine()) {
    opts.setUtmZone(FLAGS_utm_zone_num.Get(), FLAGS_utm_zone_letter.Get()[0]);
  }

  std::cout << "convert options: " << std::endl << opts << std::endl;

  SegyFile segyfile(infile);
  segyfile.open(std::ios_base::in);

  StackFile stkFile(outfile, segyfile, opts);

  segyfile.close();

  std::cout << "Num inlines: " << stkFile.getNumInlines() << std::endl;
  std::cout << "Num crosslines: " << stkFile.getNumCrosslines() << std::endl;
  std::cout << stkFile.grid() << std::endl;

  if (FLAGS_enable_crossline_opt.IsSpecifiedOnCommandLine())
    stkFile.setCrosslineAccessOptimization(FLAGS_enable_crossline_opt.Get());

  if (FLAGS_enable_depth_opt.IsSpecifiedOnCommandLine())
    stkFile.setDepthSliceAccessOptimization(FLAGS_enable_depth_opt.Get());

  return 0;
}
