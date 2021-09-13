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
ABSL_FLAG(std::string, output_stats_file, "", "Output Statistics file");
ABSL_FLAG(int, il_offset, 189, "Inline number byte offset");
ABSL_FLAG(int, xl_offset, 193, "Crossline number byte offset");
ABSL_FLAG(int, x_coord_offset, 181, "X coordinate byte offset");
ABSL_FLAG(int, y_coord_offset, 185, "Y coordinate byte offset");
ABSL_FLAG(int, utm_zone_num, 32, "UTM Zone number");
ABSL_FLAG(bool, is_2d, false, "The input file is a 2D dataset.");
ABSL_FLAG(std::string, utm_zone_letter, "U", "UTM Zone letter");
ABSL_FLAG(bool, enable_crossline_opt, false, "Enable fast crossline access");
ABSL_FLAG(bool, enable_depth_opt, false, "Enable fast depth slice access");
ABSL_FLAG(bool, verbose, false, "Enable verbose logging");
ABSL_FLAG(int, verbose_level, 0, "Verbose logging level");

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

  if (FLAGS_verbose.Get()) {
    FLAGS_stderrthreshold = 0;
  }

  if (FLAGS_verbose_level.IsSpecifiedOnCommandLine()) {
    FLAGS_v = FLAGS_verbose_level.Get();
  }

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

  SegyFile segyfile(infile);
  segyfile.open(std::ios_base::in);

  std::cout << segyfile.getTextHeader() << std::endl;
  std::cout << segyfile.getBinaryHeader() << std::endl;

  StackFile::SegyOptions opts;
  opts.setTraceHeaderOffsets(segyfile.guessTraceHeaderOffsets());

  if (FLAGS_il_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(SegyFile::Trace::Header::Attribute::INLINE_NUMBER,
                              FLAGS_il_offset.Get());
  }

  if (FLAGS_xl_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(
        SegyFile::Trace::Header::Attribute::CROSSLINE_NUMBER,
        FLAGS_xl_offset.Get());
  }

  if (FLAGS_x_coord_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(SegyFile::Trace::Header::Attribute::X_COORDINATE,
                              FLAGS_x_coord_offset.Get());
  }

  if (FLAGS_y_coord_offset.IsSpecifiedOnCommandLine()) {
    opts.setTraceHeaderOffset(SegyFile::Trace::Header::Attribute::Y_COORDINATE,
                              FLAGS_y_coord_offset.Get());
  }

  if (FLAGS_utm_zone_num.IsSpecifiedOnCommandLine() &&
      FLAGS_utm_zone_letter.IsSpecifiedOnCommandLine()) {
    opts.setUtmZone(FLAGS_utm_zone_num.Get(), FLAGS_utm_zone_letter.Get()[0]);
  }

  if (FLAGS_is_2d.Get()) {
    opts.setIs2D(FLAGS_is_2d.Get());
  }

  if (FLAGS_output_stats_file.IsSpecifiedOnCommandLine()) {
    opts.setOutputStatsFile(FLAGS_output_stats_file.Get());
  }

  std::cout << std::endl
            << "convert options: " << std::endl
            << opts << std::endl;

  StackFile stkFile(outfile, segyfile, opts);

  segyfile.close();

  std::cout << std::endl
            << "Grid: " << std::endl
            << stkFile.grid() << std::endl;

  if (FLAGS_enable_crossline_opt.Get()) {
    std::cout << "Creating crossline access optimization ..." << std::endl;
    stkFile.setCrosslineAccessOptimization(FLAGS_enable_crossline_opt.Get());
    std::cout << "Creating crossline access optimization done." << std::endl;
  }

  if (FLAGS_enable_depth_opt.Get()) {
    std::cout << "Creating depth access optimization ..." << std::endl;
    stkFile.setDepthSliceAccessOptimization(FLAGS_enable_depth_opt.Get());
    std::cout << "Creating depth access optimization done." << std::endl;
  }

  return 0;
}
