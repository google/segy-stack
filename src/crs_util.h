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

#ifndef SEGYSTACK_CRS_UTIL_H_
#define SEGYSTACK_CRS_UTIL_H_

#include <memory>
#include <ostream>
#include <set>
#include <utility>

namespace segystack {

class UTMZone {
 public:
  UTMZone(int zone_num = 32, char zone_char = 'U');  // default UTM Zone.

  std::pair<int, char> value() const;
  void setValue(int zone_num, char zone_char);

  int number() const;
  void setNumber(int zone_num);

  char letter() const;
  void setLetter(char zone_char);

 private:
  int zone_num_;
  char zone_char_;
};

struct GeographicCoordinates {
  float latitude;
  float longitude;
};

class UTMZoneConverter {
 public:
  UTMZoneConverter(const UTMZone& utm_zone);

  GeographicCoordinates getGeographicCoordinates(float easting,
                                                 float northing) const;

  ~UTMZoneConverter();
 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

std::ostream& operator<<(std::ostream& os, const UTMZone& utm);

}  // namespace segystack

#endif
