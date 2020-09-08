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

#ifndef SEGYSTACK_SEGY_H_
#define SEGYSTACK_SEGY_H_

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "float_util.h"
#include "mmap_file.h"

#ifdef __APPLE__
#include <machine/endian.h>
#else
#include <endian.h>
#endif

namespace segystack {

// Reads and writes SEGY files.
// http://seg.org/Portals/0/SEG/News%20and%20Resources/Technical%20Standards/seg_y_rev1.pdf
class SegyFile {
 public:
  static const size_t TEXT_HEADER_SIZE = 3200;
  static const size_t BINARY_HEADER_SIZE = 400;
  static const size_t TRACE_HEADER_SIZE = 240;

  // The first 3200-byte, Textual File Header record contains 40 lines of
  // textual information, providing a human-readable description of the seismic
  // data in the SEG-Y file.
  struct TextHeader {
    char data[TEXT_HEADER_SIZE];
    explicit operator std::string() const { return this->toString(); }
    std::string toString() const;
    void set_value(const std::string& value);
  };

// Returns value of specified type at 1 numbered byte |offset|.
#define ADD_VALUE_AT_OFFSET_GETTER_METHOD                          \
  template <typename T>                                            \
  T getValueAtOffset(int offset) const {                           \
    offset--;                                                      \
    const char* addr = data;                                       \
    const T* val_addr = reinterpret_cast<const T*>(addr + offset); \
    T value = *val_addr;                                           \
    value = fix_endianness_if_needed(value);                       \
    return value;                                                  \
  }

  // The 400-byte Binary File Header record contains binary values that affect
  // the whole SEG Y file.
  struct BinaryHeader {
    enum Attribute {
      JOB_IDENT_NUM,
      LINE_NUM,
      REEL_NUM,
      NUM_DATA_TRC_ENSEMBLE,
      NUM_AUX_TRC_ENSEMBLE,
      SAMP_INT,
      SAMP_INT_FIELD,
      NUM_SAMP_PER_TRC,
      NUM_SAMP_PER_TRC_FIELD,
      SAMPLE_FORMAT_CODE
    };

    char data[BINARY_HEADER_SIZE];

    ADD_VALUE_AT_OFFSET_GETTER_METHOD
    void print(std::ostream& os = std::cout) const;
  };

  class Trace {
   public:
    // The SEG Y trace header contains trace attributes.
    class Header {
     public:
      enum class Attribute {
        INLINE_NUMBER,
        CROSSLINE_NUMBER,
        ENSEMBLE_NUMBER,
        SHOTPOINT_NUMBER,
        X_COORDINATE,
        Y_COORDINATE
      };

      char data[TRACE_HEADER_SIZE];

      ADD_VALUE_AT_OFFSET_GETTER_METHOD

      // Applies scalar corrections to a integer encoded coordinate value and
      // returns a floating point value.
      double getCoordinateValue(int offset) const;

      void print(std::ostream& os = std::cout) const;
    };

    const Header& header() const { return header_; }
    Header& header() { return header_; }

    const std::vector<float>& samples() const { return samples_; }
    std::vector<float>& samples() { return samples_; }

   private:
    Header header_;
    std::vector<float> samples_;
  };

  explicit SegyFile(const std::string& filename);

  void open(std::ios_base::openmode mode = std::ios_base::in);

  std::ios_base::openmode open_mode() const { return mode_; }

  bool is_open() const { return file_->is_open(); }

  const std::string& name() const { return file_->name(); }

  void close();

  std::map<Trace::Header::Attribute, int> guessTraceHeaderOffsets() const;

  const TextHeader& getTextHeader() const;
  void setTextHeader(const TextHeader& header);

  const BinaryHeader& getBinaryHeader() const;
  void setBinaryHeader(const BinaryHeader& header);

  // advances both the trace header and sample pointers
  // to the offset from the beginning of the file.
  void seek(std::uint64_t offset) const;

  // reads a Trace but does not advance the internal trace
  // pointer to the next trace header. Must seek to the correct offset
  // before calling this method
  bool read(Trace& trace) const;

 protected:
  void checkFileOpened() const;
  void checkFileOpenedForReading() const;

  std::unique_ptr<MmapFile> file_;
  std::ios_base::openmode mode_;
  TextHeader* text_header_;
  BinaryHeader* binary_header_;
  std::uint16_t num_samples_per_trc_;
  std::int16_t num_ext_hdrs_;
  char* first_hdr_ptr_;
  char* hdr_ptr_;
  char* trc_ptr_;
  std::uint64_t cur_offset_;
};

}  // namespace segystack

std::ostream& operator<<(std::ostream& os,
                         const segystack::SegyFile::TextHeader& hdr);

std::ostream& operator<<(std::ostream& os,
                         const segystack::SegyFile::BinaryHeader& hdr);

std::ostream& operator<<(std::ostream& os,
                         const segystack::SegyFile::Trace::Header& hdr);

#endif  // SEGYSTACK_SEGY_H_
