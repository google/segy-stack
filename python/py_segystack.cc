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

#include <glog/logging.h>

#include <sstream>

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "segy_file.h"
#include "stack_file.h"
#include "test/test_base.h"

namespace py = pybind11;
using namespace py::literals;
using namespace segystack;

void init_common(py::module* m) {
  py::class_<UTMZone> utm(*m, "UTMZone");
  utm.def(py::init<>());
  utm.def(py::init<int, char>());
  utm.def_property("number", &UTMZone::number, &UTMZone::setNumber);
  utm.def_property(
      "letter", &UTMZone::letter, [](UTMZone& self, const std::string& val) {
        if (val.size() != 1)
          throw py::type_error("Zone must consist of a single character!");
        self.setLetter(val[0]);
      });
  utm.def("__repr__", [](const UTMZone& val) {
    std::ostringstream ostr;
    ostr << val;
    return ostr.str();
  });
}

void init_segy_file(py::module* m) {
  py::class_<SegyFile> sgy(*m, "SegyFile");

  sgy.def(py::init<const std::string&>());
  sgy.def("open", [](SegyFile& self, const std::string& mode) {
    std::ios_base::openmode open_mode = std::ios_base::in;
    if (mode.find('r') != std::string::npos)
      open_mode |= std::ios_base::in;
    if (mode.find('w') != std::string::npos)
      open_mode |= std::ios_base::out;
    self.open(open_mode);
  });
  sgy.def("name", &SegyFile::name);
  sgy.def("close", &SegyFile::close);
  sgy.def("guess_trace_header_offsets", &SegyFile::guessTraceHeaderOffsets);

  py::class_<SegyFile::Trace> trace(sgy, "Trace");
  trace.def(py::init<>());
  trace.def("header",
            [](const SegyFile::Trace& self) { return self.header(); });

  py::class_<SegyFile::Trace::Header> hdr(trace, "Header");
  hdr.def(py::init<>());

  py::enum_<SegyFile::Trace::Header::Attribute>(hdr, "Attribute")
      .value("INLINE_NUMBER", SegyFile::Trace::Header::Attribute::INLINE_NUMBER)
      .value("CROSSLINE_NUMBER",
             SegyFile::Trace::Header::Attribute::CROSSLINE_NUMBER)
      .value("X_COORDINATE", SegyFile::Trace::Header::Attribute::X_COORDINATE)
      .value("Y_COORDINATE", SegyFile::Trace::Header::Attribute::Y_COORDINATE)
      .value("SHOTPOINT_NUMBER", SegyFile::Trace::Header::Attribute::SHOTPOINT_NUMBER)
      .value("ENSEMBLE_NUMBER", SegyFile::Trace::Header::Attribute::ENSEMBLE_NUMBER)
      .export_values();
}

void init_stack_file(py::module* m) {
  py::class_<StackFile> sf(*m, "StackFile");

  sf.def(py::init<const std::string&, const SegyFile&,
                  const StackFile::SegyOptions&>());
  sf.def(py::init<const std::string&, const SegyFile&>());
  sf.def(py::init<const std::string&>());
  sf.def("name", &StackFile::name);
  sf.def("grid", &StackFile::grid, py::return_value_policy::reference);
  sf.def("read_inline",
         [](const StackFile& self, int il,
            float fill_value) -> py::array_t<float> {
           py::array::ShapeContainer shape(
               {self.grid().numCrosslines(), self.grid().numSamples()});
           py::array_t<float> data(shape);
           float* buffer = data.mutable_data();
           self.readInline(il, buffer, data.size(), fill_value);
           return data;
         });
  sf.def("read_crossline",
         [](const StackFile& self, int xl,
            float fill_value) -> py::array_t<float> {
           py::array::ShapeContainer shape(
               {self.grid().numInlines(), self.grid().numSamples()});
           py::array_t<float> data(shape);
           float* buffer = data.mutable_data();
           self.readCrossline(xl, buffer, data.size(), fill_value);
           return data;
         });
  sf.def("read_depth_slice",
         [](const StackFile& self, int sample_index,
            float fill_value) -> py::array_t<float> {
           py::array::ShapeContainer shape(
               {self.grid().numInlines(), self.grid().numCrosslines()});
           py::array_t<float> data(shape);
           float* buffer = data.mutable_data();
           self.readDepthSlice(sample_index, buffer, data.size(), fill_value);
           return data;
         });
  sf.def("set_crossline_access_opt",
         &StackFile::setCrosslineAccessOptimization);
  sf.def("set_depth_slice_access_opt",
         &StackFile::setDepthSliceAccessOptimization);

  py::class_<StackFile::Grid> grid(sf, "Grid");
  grid.def(py::init<>());
  grid.def("utm_zone", &StackFile::Grid::utmZone);
  grid.def("bounding_box", &StackFile::Grid::boundingBox);
  grid.def("is_bin_active", &StackFile::Grid::isBinActive);
  grid.def("get_coordinate", &StackFile::Grid::getCoordinate);
  grid.def_property("inline_min", &StackFile::Grid::inlineMin,
                    &StackFile::Grid::setInlineMin);
  grid.def_property("inline_max", &StackFile::Grid::inlineMax,
                    &StackFile::Grid::setInlineMax);
  grid.def_property("inline_increment", &StackFile::Grid::inlineIncrement,
                    &StackFile::Grid::setInlineIncrement);
  grid.def_property("inline_spacing", &StackFile::Grid::inlineSpacing,
                    &StackFile::Grid::setInlineSpacing);
  grid.def_property("crossline_min", &StackFile::Grid::crosslineMin,
                    &StackFile::Grid::setCrosslineMin);
  grid.def_property("crossline_max", &StackFile::Grid::crosslineMax,
                    &StackFile::Grid::setCrosslineMax);
  grid.def_property("crossline_increment", &StackFile::Grid::crosslineIncrement,
                    &StackFile::Grid::setCrosslineIncrement);
  grid.def_property("crossline_spacing", &StackFile::Grid::crosslineSpacing,
                    &StackFile::Grid::setCrosslineSpacing);
  grid.def_property("num_inlines", &StackFile::Grid::numInlines,
                    &StackFile::Grid::setNumInlines);
  grid.def_property("num_crosslines", &StackFile::Grid::numCrosslines,
                    &StackFile::Grid::setNumCrosslines);
  grid.def_property("sampling_interval", &StackFile::Grid::samplingInterval,
                    &StackFile::Grid::setSamplingInterval);
  grid.def_property("num_samples", &StackFile::Grid::numSamples,
                    &StackFile::Grid::setNumSamples);
  grid.def_property("units", &StackFile::Grid::units,
                    &StackFile::Grid::setUnits);
  grid.def_property(
      "inline_numbers",
      [](const StackFile::Grid& self) {
        py::tuple ils(self.numInlines());
        for (int il = self.inlineMin(), idx = 0; il <= self.inlineMax();
             il += self.inlineIncrement(), ++idx) {
          ils[idx] = il;
        }
        return ils;
      },
      [](const StackFile::Grid&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  grid.def_property(
      "crossline_numbers",
      [](const StackFile::Grid& self) {
        py::tuple xls(self.numCrosslines());
        for (int xl = self.crosslineMin(), idx = 0; xl <= self.crosslineMax();
             xl += self.crosslineIncrement(), ++idx) {
          xls[idx] = xl;
        }
        return xls;
      },
      [](const StackFile::Grid&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  grid.def("__repr__", [](const StackFile::Grid& grid) {
    std::ostringstream ostr;
    ostr << grid;
    return ostr.str();
  });

  py::enum_<StackFile::Grid::Units>(grid, "Units")
      .value("METERS", StackFile::Grid::METERS)
      .value("FEET", StackFile::Grid::FEET)
      .export_values();

  py::class_<StackFile::Grid::Coordinate> coord(grid, "Coordinate");
  coord.def(py::init<>());
  coord.def_property(
      "x", [](const StackFile::Grid::Coordinate& self) { return self.x; },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "y", [](const StackFile::Grid::Coordinate& self) { return self.y; },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "inline_num",
      [](const StackFile::Grid::Coordinate& self) { return self.inline_num; },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "crossline_num",
      [](const StackFile::Grid::Coordinate& self) {
        return self.crossline_num;
      },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "ensemble_num",
      [](const StackFile::Grid::Coordinate& self) { return self.ensemble_num; },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "shotpoint_num",
      [](const StackFile::Grid::Coordinate& self) {
        return self.shotpoint_num;
      },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "latitude",
      [](const StackFile::Grid::Coordinate& self) { return self.lat; },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def_property(
      "longitude",
      [](const StackFile::Grid::Coordinate& self) { return self.lon; },
      [](StackFile::Grid::Coordinate&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  coord.def("__repr__", [](const StackFile::Grid::Coordinate& val) {
    std::ostringstream ostr;
    ostr << val;
    return ostr.str();
  });

  py::class_<StackFile::Grid::BoundingBox> bbox(grid, "BoundingBox");
  bbox.def(py::init<>());
  bbox.def_property(
      "c1", [](const StackFile::Grid::BoundingBox& self) { return self.c1; },
      [](StackFile::Grid::BoundingBox&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  bbox.def_property(
      "c2", [](const StackFile::Grid::BoundingBox& self) { return self.c2; },
      [](StackFile::Grid::BoundingBox&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  bbox.def_property(
      "c3", [](const StackFile::Grid::BoundingBox& self) { return self.c3; },
      [](StackFile::Grid::BoundingBox&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  bbox.def_property(
      "c4", [](const StackFile::Grid::BoundingBox& self) { return self.c4; },
      [](StackFile::Grid::BoundingBox&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  bbox.def("__repr__", [](const StackFile::Grid::BoundingBox& val) {
    std::ostringstream ostr;
    ostr << val;
    return ostr.str();
  });

  py::class_<StackFile::SegyOptions> opts(sf, "SegyOptions");
  opts.def(py::init<>());
  opts.def("set_utm_zone",
           [](StackFile::SegyOptions& self, int number, char v) {
             try {
               self.setUtmZone(number, v);
             } catch (const std::runtime_error& e) {
               throw py::type_error(e.what());
             }
           });
  opts.def("utm_zone", &StackFile::SegyOptions::getUtmZone);
  opts.def("set_trace_header_offset",
           &StackFile::SegyOptions::setTraceHeaderOffset);
  opts.def("set_trace_header_offsets",
           &StackFile::SegyOptions::setTraceHeaderOffsets);
  opts.def("trace_header_offset",
           &StackFile::SegyOptions::getTraceHeaderOffset);
  opts.def("set_is_2D", &StackFile::SegyOptions::setIs2D);
  opts.def("is_2D", &StackFile::SegyOptions::is2D);
  opts.def("set_output_stats_file", &StackFile::SegyOptions::setOutputStatsFile);
  opts.def("output_stats_file", &StackFile::SegyOptions::getOutputStatsFile);
  opts.def("__repr__", [](const StackFile::SegyOptions& opts) {
    std::ostringstream ostr;
    ostr << opts;
    return ostr.str();
  });
}

PYBIND11_MODULE(segystack, m) {
  google::InitGoogleLogging("segystack");

  m.doc() = "segystack python interface";
  init_common(&m);
  init_segy_file(&m);
  init_stack_file(&m);

  py::module test_mod = m.def_submodule("test", "Testing utility methods");
  test_mod.def(
      "create_test_segy", &segystack::test::create_test_segy,
      "Helper function to create test SEG-Y files.", py::arg("outfile"),
      py::arg("num_samples"), py::arg("sampling_interval"), py::arg("num_il"),
      py::arg("il_increment"), py::arg("num_xl"), py::arg("xl_increment"),
      py::arg("x_origin"), py::arg("y_origin"), py::arg("il_spacing"),
      py::arg("xl_spacing"), py::arg("opts"));
}
