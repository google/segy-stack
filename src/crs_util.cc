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

#include "crs_util.h"

#include "config.h"
#include "logging.h"
#include "proj.h"

#include <sstream>
#include <vector>

namespace segystack {

std::ostream& operator<<(std::ostream& os, const UTMZone& utm) {
  os << utm.value().first << " " << utm.value().second;
  return os;
}

UTMZone::UTMZone(int zone_num, char zone_char) {
  setValue(zone_num, zone_char);
}

std::pair<int, char> UTMZone::value() const {
  return std::make_pair(zone_num_, zone_char_);
}

int UTMZone::number() const {
  return zone_num_;
}

void UTMZone::setNumber(int zone_num) {
  setValue(zone_num, zone_char_);
}

char UTMZone::letter() const {
  return zone_char_;
}

void UTMZone::setLetter(char zone_char) {
  setValue(zone_num_, zone_char);
}

void UTMZone::setValue(int zone_num, char zone_char) {
  if (zone_num < 1 || zone_num > 60)
    throw std::runtime_error("UTM zone number should be in range [1, 60]");

  zone_char = std::toupper(zone_char);
  std::set<char> zones;
  const std::set<char> forbidden = {'A', 'B', 'Y', 'Z', 'I', 'O'};
  for (char c = 'A'; c <= 'Z'; c++) {
    if (forbidden.find(c) == forbidden.end())
      zones.insert(c);
  }

  if (zones.find(zone_char) == zones.end()) {
    std::ostringstream ostr;
    for (char c : zones)
      ostr << c << ", ";
    throw std::runtime_error("UTM zone letter should be one of " + ostr.str());
  }

  zone_num_ = zone_num;
  zone_char_ = zone_char;
}

class UTMZoneConverter::Impl {
 public:
  Impl(const UTMZone& utm) {
    std::vector<std::string> path_strs;
    std::vector<const char*> path_cstrs;

    std::string proj_data_paths(SEGYSTACK_PROJ_DIR);

    size_t start = 0;
    size_t end = proj_data_paths.find(':');
    while (end != std::string::npos) {
      std::string proj_dir = proj_data_paths.substr(start, end - start);
      LOGFN_VAR(proj_dir);
      path_strs.push_back(proj_dir);
      path_cstrs.push_back(path_strs.back().c_str());
      start = end + 1;
      end = proj_data_paths.find(':', start);
    }
    std::string proj_dir = proj_data_paths.substr(start);
    LOGFN_VAR(proj_dir);
    path_strs.push_back(proj_dir);
    path_cstrs.push_back(path_strs.back().c_str());

    proj_context_set_search_paths(PJ_DEFAULT_CTX, path_cstrs.size(),
                                  &path_cstrs[0]);

    std::ostringstream proj_str;
    proj_str << "+proj=utm +zone=" << utm.number();
    if (utm.letter() < 'N') {
      proj_str << " +south";
    }
    proj_str << " +datum=WGS84";
    LOGFN_VAR(proj_str.str());

    PJ* pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326",
                                    proj_str.str().c_str(), NULL);

    if (!pj) {
      LOG(WARNING) << "Could not create PROJ context! Latitude & Longitude "
                      "info will be unavailable";
    }
    proj_context_ = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
    proj_destroy(pj);

    if (!proj_context_) {
      LOG(WARNING) << "Could not create PROJ context! Latitude & Longitude "
                      "info will be unavailable";
    }
  }

  ~Impl() { proj_destroy(proj_context_); }

 private:
  friend class UTMZoneConverter;
  PJ* proj_context_;
};

UTMZoneConverter::UTMZoneConverter(const UTMZone& utm_zone)
    : impl_(new Impl(utm_zone)) {}

UTMZoneConverter::~UTMZoneConverter() = default;

GeographicCoordinates UTMZoneConverter::getGeographicCoordinates(
    float easting,
    float northing) const {
  if (impl_->proj_context_ == nullptr)
    return GeographicCoordinates();

  PJ_COORD a;
  a.xy.x = easting;
  a.xy.y = northing;
  PJ_COORD b = proj_trans(impl_->proj_context_, PJ_INV, a);

  GeographicCoordinates coord;
  coord.longitude = b.lp.lam;
  coord.latitude = b.lp.phi;

  return coord;
}

}  // namespace segystack
