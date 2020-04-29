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
#include "stack_file.h"
#include "test/test_base.h"

namespace py = pybind11;
using namespace py::literals;
using namespace segystack;

void init_types(py::module* m) {
  py::class_<Grid> grid(*m, "Grid");
  grid.def(py::init<>());
  grid.def_property("inline_min", &Grid::inline_min, &Grid::set_inline_min);
  grid.def_property("inline_max", &Grid::inline_max, &Grid::set_inline_max);
  grid.def_property("inline_increment", &Grid::inline_increment,
                    &Grid::set_inline_increment);
  grid.def_property("inline_spacing", &Grid::inline_spacing,
                    &Grid::set_inline_spacing);
  grid.def_property("crossline_min", &Grid::crossline_min,
                    &Grid::set_crossline_min);
  grid.def_property("crossline_max", &Grid::crossline_max,
                    &Grid::set_crossline_max);
  grid.def_property("crossline_increment", &Grid::crossline_increment,
                    &Grid::set_crossline_increment);
  grid.def_property("crossline_spacing", &Grid::crossline_spacing,
                    &Grid::set_crossline_spacing);
  grid.def_property("sampling_interval", &Grid::sampling_interval,
                    &Grid::set_sampling_interval);
  grid.def_property("num_samples", &Grid::num_samples, &Grid::set_num_samples);
  grid.def_property("units", &Grid::units, &Grid::set_units);
  grid.def_property("num_active_cells", &Grid::num_active_cells,
                    &Grid::set_num_active_cells);
  grid.def("__repr__", [](const Grid& grid) {
    std::ostringstream ostr;
    ostr << grid;
    return ostr.str();
  });

  py::enum_<Grid::Units>(grid, "Units")
      .value("METERS", Grid::METERS)
      .value("FEET", Grid::FEET)
      .export_values();

  py::class_<Grid::Cell> cell(grid, "Cell");
  cell.def(py::init<>());
  cell.def_property("x_coordinate", &Grid::Cell::x_coordinate,
                    &Grid::Cell::set_x_coordinate);
  cell.def_property("y_coordinate", &Grid::Cell::y_coordinate,
                    &Grid::Cell::set_y_coordinate);
  cell.def_property("inline_number", &Grid::Cell::inline_number,
                    &Grid::Cell::set_inline_number);
  cell.def_property("crossline_number", &Grid::Cell::crossline_number,
                    &Grid::Cell::set_crossline_number);
  cell.def("__repr__", [](const Grid::Cell& cell) {
    std::ostringstream ostr;
    ostr << cell;
    return ostr.str();
  });

  py::class_<UTMZone> utm(*m, "UTMZone");
  utm.def(py::init<>());
  utm.def_property("number", &UTMZone::number, &UTMZone::set_number);
  utm.def_property(
      "letter", &UTMZone::letter,
      [](UTMZone& self, const std::string& val) { self.set_letter(val); });
  utm.def("__repr__", [](const UTMZone& val) {
    std::ostringstream ostr;
    ostr << val;
    return ostr.str();
  });
}

void init_stack_file(py::module* m) {
  py::class_<StackFile> sf(*m, "StackFile");

  sf.def(py::init<const std::string&, const std::string&,
                  const StackFile::SegyOptions&>());
  sf.def(py::init<const std::string&>());
  sf.def("grid", &StackFile::grid, py::return_value_policy::reference);
  sf.def("num_inlines", &StackFile::getNumInlines);
  sf.def("num_crosslines", &StackFile::getNumCrosslines);
  sf.def("utm_zone", &StackFile::getUtmZone);
  sf.def_property(
      "inline_numbers",
      [](const StackFile& self) {
        py::tuple ils(self.getNumInlines());
        for (int il = self.grid().inline_min(), idx = 0;
             il <= self.grid().inline_max();
             il += self.grid().inline_increment(), ++idx) {
          ils[idx] = il;
        }
        return ils;
      },
      [](const StackFile&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  sf.def_property(
      "crossline_numbers",
      [](const StackFile& self) {
        py::tuple xls(self.getNumCrosslines());
        for (int xl = self.grid().crossline_min(), idx = 0;
             xl <= self.grid().crossline_max();
             xl += self.grid().crossline_increment(), ++idx) {
          xls[idx] = xl;
        }
        return xls;
      },
      [](const StackFile&, const py::object&) {
        throw py::type_error("Read only attribute");
      });
  sf.def("read_inline",
         [](const StackFile& self, int il,
            float fill_value) -> py::array_t<float> {
           py::array::ShapeContainer shape(
               {self.getNumCrosslines(), self.grid().num_samples()});
           py::array_t<float> data(shape);
           float* buffer = data.mutable_data();
           self.readInline(il, buffer, data.size(), fill_value);
           return data;
         });
  sf.def("read_crossline",
         [](const StackFile& self, int xl,
            float fill_value) -> py::array_t<float> {
           py::array::ShapeContainer shape(
               {self.getNumInlines(), self.grid().num_samples()});
           py::array_t<float> data(shape);
           float* buffer = data.mutable_data();
           self.readCrossline(xl, buffer, data.size(), fill_value);
           return data;
         });
  sf.def("read_depth_slice",
         [](const StackFile& self, int sample_index,
            float fill_value) -> py::array_t<float> {
           py::array::ShapeContainer shape(
               {self.getNumInlines(), self.getNumCrosslines()});
           py::array_t<float> data(shape);
           float* buffer = data.mutable_data();
           self.readDepthSlice(sample_index, buffer, data.size(), fill_value);
           return data;
         });
  sf.def("set_crossline_access_opt",
         &StackFile::setCrosslineAccessOptimization);
  sf.def("set_depth_slice_access_opt",
         &StackFile::setDepthSliceAccessOptimization);

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
  opts.def("trace_header_offset",
           &StackFile::SegyOptions::getTraceHeaderOffset);
  opts.def("__repr__", [](const StackFile::SegyOptions& opts) {
    std::ostringstream ostr;
    ostr << opts;
    return ostr.str();
  });

  py::enum_<StackFile::SegyOptions::TraceHeaderAttribute>(
      opts, "TraceHeaderAttribute")
      .value("INLINE_NUMBER", StackFile::SegyOptions::INLINE_NUMBER)
      .value("CROSSLINE_NUMBER", StackFile::SegyOptions::CROSSLINE_NUMBER)
      .value("X_COORDINATE", StackFile::SegyOptions::X_COORDINATE)
      .value("Y_COORDINATE", StackFile::SegyOptions::Y_COORDINATE)
      .export_values();
}

PYBIND11_MODULE(segystack, m) {
  google::InitGoogleLogging("segystack");

  m.doc() = "segystack python interface";
  init_types(&m);
  init_stack_file(&m);

  py::module test_mod = m.def_submodule("test", "Testing utility methods");
  test_mod.def("create_test_segy", &segystack::test::create_test_segy);
}
