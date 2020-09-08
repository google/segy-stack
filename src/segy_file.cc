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

#include "segy_file.h"

#include <algorithm>
#include <cstring>
#include <map>
#include <sstream>
#include <utility>

#include "logging.h"

std::ostream& operator<<(std::ostream& os,
                         const segystack::SegyFile::TextHeader& hdr) {
  os << hdr.toString();
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const segystack::SegyFile::BinaryHeader& hdr) {
  hdr.print(os);
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const segystack::SegyFile::Trace::Header& hdr) {
  hdr.print(os);
  return os;
}

namespace segystack {

namespace {

constexpr int kSegyNumSamplesPerTrcOffset = 21;
constexpr int kSegyNumExtendedHdsOffset = 305;
constexpr int kSegyIEEEDataFormatCode = 5;
constexpr int kSegyIBMFloatDataFormatCode = 1;
constexpr int kSegyDataSampleFormat = 25;
constexpr int kSegyTextHeaderLineWidth = 80;
constexpr int kSegyBinHdrSamplingIntOffset = 17;
constexpr int kSegyBinHdrSegyMajorOffset = 301;
constexpr int kSegyBinHdrSegyMinorOffset = 302;
constexpr int kSegyHeaderEnsembleNumOffset = 21;
constexpr int kSegyHeaderShotpointNumOffset = 197;

const std::vector<int> kSegyHeaderXCoordCandidateOffsets = {
    181,  // X coordinate of ensemble (CDP) position of this trace.
    73,   // Source coordinate - X.
    201,  // Non-standard but often used.
};

const std::vector<int> kSegyHeaderYCoordCandidateOffsets = {
    185,  // Y coordinate of ensemble (CDP) position of this trace.
    77,   // Source coordinate - Y.
    205,  // Non-standard but often used.
};

const std::vector<int> kSegyHeaderInlineCandidateOffsets = {
    189,  // For 3-D poststack data, this field should be used for the in-line
          // number.
    213,  // Non-standard but often used.
    9,    // Original field record number.
};

const std::vector<int> kSegyHeaderCrosslineCandidateOffsets = {
    193,  // For 3-D poststack data, this field should be used for the
          // cross-line number.
    217,  // Non-standard but often used.
    21,   // Ensemble number (i.e. CDP, CMP, CRP, etc).
};

// https://www.ibm.com/support/knowledgecenter/SS2MB5_14.1.0/com.ibm.xlf141.bg.doc/language_ref/asciit.html
const std::map<unsigned char, char> kEBCDICtoASCIImap = {
    {75, '.'},  {76, '<'},   {77, '('},  {78, '+'},  {79, '|'},  {80, '&'},
    {90, '!'},  {91, '$'},   {92, '*'},  {93, ')'},  {94, ';'},  {96, '-'},
    {97, '/'},  {106, '|'},  {107, ','}, {108, '%'}, {109, '_'}, {110, '>'},
    {111, '?'}, {121, '`'},  {122, ':'}, {123, '#'}, {124, '@'}, {125, '\''},
    {126, '='}, {127, '"'},  {129, 'a'}, {130, 'b'}, {131, 'c'}, {132, 'd'},
    {133, 'e'}, {134, 'f'},  {135, 'g'}, {136, 'h'}, {137, 'i'}, {145, 'j'},
    {146, 'k'}, {147, 'l'},  {148, 'm'}, {149, 'n'}, {150, 'o'}, {151, 'p'},
    {152, 'q'}, {153, 'r'},  {161, '~'}, {162, 's'}, {163, 't'}, {164, 'u'},
    {165, 'v'}, {166, 'w'},  {167, 'x'}, {168, 'y'}, {169, 'z'}, {192, '{'},
    {193, 'A'}, {194, 'B'},  {195, 'C'}, {196, 'D'}, {197, 'E'}, {198, 'F'},
    {199, 'G'}, {200, 'H'},  {201, 'I'}, {208, '}'}, {209, 'J'}, {210, 'K'},
    {211, 'L'}, {212, 'M'},  {213, 'N'}, {214, 'O'}, {215, 'P'}, {216, 'Q'},
    {217, 'R'}, {224, '\\'}, {226, 'S'}, {227, 'T'}, {228, 'U'}, {229, 'V'},
    {230, 'W'}, {231, 'X'},  {232, 'Y'}, {233, 'Z'}, {240, '0'}, {241, '1'},
    {242, '2'}, {243, '3'},  {244, '4'}, {245, '5'}, {246, '6'}, {247, '7'},
    {248, '8'}, {249, '9'}};

char getASCIIForEBCDIC(char c) {
  if (kEBCDICtoASCIImap.find(c) != kEBCDICtoASCIImap.end())
    return kEBCDICtoASCIImap.at(c);
  return ' ';
}

bool isTextInEBCDICFormat(const char* text, size_t length) {
  int alnumASCII = 0;
  for (size_t i = 0; i < length; i++) {
    if (std::isalnum(text[i]))
      alnumASCII++;
  }

  int alnumEBCDIC = 0;
  for (size_t i = 0; i < length; i++) {
    if (std::isalnum(getASCIIForEBCDIC(text[i])))
      alnumEBCDIC++;
  }

  if (alnumASCII > alnumEBCDIC)
    return false;
  return true;
}
}  // namespace

const size_t SegyFile::TEXT_HEADER_SIZE;
const size_t SegyFile::BINARY_HEADER_SIZE;
const size_t SegyFile::TRACE_HEADER_SIZE;

SegyFile::SegyFile(const std::string& filename)
    : text_header_(nullptr),
      binary_header_(nullptr),
      num_samples_per_trc_(0),
      first_hdr_ptr_(nullptr),
      hdr_ptr_(nullptr),
      trc_ptr_(nullptr),
      cur_offset_(0) {
  file_ = MmapFile::Create(filename);
}

std::map<SegyFile::Trace::Header::Attribute, int>
SegyFile::guessTraceHeaderOffsets() const {
  if (!is_open()) {
    throw std::runtime_error("File " + name() + " not opened for reading!");
  }

  std::uint64_t prev_offset = cur_offset_;

  Trace trace1, trace2;
  seek(0);
  read(trace1);
  seek(1);
  read(trace2);

  // restore back to where we were.
  seek(prev_offset);

  const Trace::Header& header1 = trace1.header();
  const Trace::Header& header2 = trace2.header();

  LOGFN_VAR(header1);
  LOGFN_VAR(header2);

  CHECK_EQ(kSegyHeaderXCoordCandidateOffsets.size(),
           kSegyHeaderYCoordCandidateOffsets.size());

  std::map<Trace::Header::Attribute, int> offsets;
  offsets[Trace::Header::Attribute::ENSEMBLE_NUMBER] =
      kSegyHeaderEnsembleNumOffset;
  offsets[Trace::Header::Attribute::SHOTPOINT_NUMBER] =
      kSegyHeaderShotpointNumOffset;

  for (size_t i = 0; i < kSegyHeaderXCoordCandidateOffsets.size(); i++) {
    int x_offset = kSegyHeaderXCoordCandidateOffsets[i];
    int y_offset = kSegyHeaderYCoordCandidateOffsets[i];

    float x_coord1 = header1.getCoordinateValue(x_offset);
    float y_coord1 = header1.getCoordinateValue(y_offset);
    float x_coord2 = header2.getCoordinateValue(x_offset);
    float y_coord2 = header2.getCoordinateValue(y_offset);

    if ((x_coord1 != x_coord2) || (y_coord1 != y_coord2)) {
      offsets[Trace::Header::Attribute::X_COORDINATE] = x_offset;
      offsets[Trace::Header::Attribute::Y_COORDINATE] = y_offset;
      break;
    }
  }

  CHECK_EQ(kSegyHeaderInlineCandidateOffsets.size(),
           kSegyHeaderCrosslineCandidateOffsets.size());

  for (size_t i = 0; i < kSegyHeaderInlineCandidateOffsets.size(); i++) {
    int il_offset = kSegyHeaderInlineCandidateOffsets[i];
    int xl_offset = kSegyHeaderCrosslineCandidateOffsets[i];

    int32_t il1 = header1.getValueAtOffset<int32_t>(il_offset);
    int32_t xl1 = header1.getValueAtOffset<int32_t>(xl_offset);
    int32_t il2 = header2.getValueAtOffset<int32_t>(il_offset);
    int32_t xl2 = header2.getValueAtOffset<int32_t>(xl_offset);

    if ((il1 != il2) || (xl1 != xl2)) {
      offsets[Trace::Header::Attribute::INLINE_NUMBER] = il_offset;
      offsets[Trace::Header::Attribute::CROSSLINE_NUMBER] = xl_offset;
      break;
    }
  }

  auto check_offset_exists = [&](Trace::Header::Attribute attr,
                                 const std::string& attr_name) {
    if (offsets.find(attr) == offsets.end()) {
      LOG(WARNING) << "Warning: Could not guess the location of " << attr_name
                   << " in the trace header!" << std::endl;
    }
  };

  check_offset_exists(Trace::Header::Attribute::INLINE_NUMBER, "inline number");
  check_offset_exists(Trace::Header::Attribute::CROSSLINE_NUMBER,
                      "crossline number");
  check_offset_exists(Trace::Header::Attribute::X_COORDINATE, "X coordinate");
  check_offset_exists(Trace::Header::Attribute::Y_COORDINATE, "Y coordinate");

  return offsets;
}

void SegyFile::open(std::ios_base::openmode mode) {
  if (file_->is_open()) {
    close();
  }

  mode_ = mode;
  file_->open(mode);
  file_->map();

  char* start_addr = file_->char_addr();
  text_header_ = reinterpret_cast<TextHeader*>(start_addr);

  LOGFN_VAR((*text_header_));

  start_addr += sizeof(TextHeader);
  binary_header_ = reinterpret_cast<BinaryHeader*>(start_addr);

  num_samples_per_trc_ =
      binary_header_->getValueAtOffset<uint16_t>(kSegyNumSamplesPerTrcOffset);
  num_ext_hdrs_ =
      binary_header_->getValueAtOffset<int16_t>(kSegyNumExtendedHdsOffset);

  LOGFN_VAR((*binary_header_));
  LOGFN_VAR(num_samples_per_trc_);
  LOGFN_VAR(num_ext_hdrs_);

  if (num_ext_hdrs_ < 0) {
    throw std::runtime_error(
        "SEGY: Cannot handle variable number of extended headers yet!");
  }

  first_hdr_ptr_ = file_->char_addr() + sizeof(TextHeader) +
                   sizeof(BinaryHeader) + num_ext_hdrs_ * sizeof(TextHeader);
  seek(0);
}

double SegyFile::Trace::Header::getCoordinateValue(int offset) const {
  double coord = double(getValueAtOffset<int32_t>(offset));
  int16_t scalar_apply_coords = getValueAtOffset<int16_t>(71);
  if (scalar_apply_coords != 0) {
    if (scalar_apply_coords > 0) {
      coord *= double(scalar_apply_coords);
    } else {
      coord /= double(std::abs(scalar_apply_coords));
    }
  }

  return coord;
}

void SegyFile::close() {
  text_header_ = nullptr;
  binary_header_ = nullptr;
  file_->unmap();
  file_->close();
}

void SegyFile::seek(std::uint64_t offset) const {
  checkFileOpened();
  SegyFile* self = const_cast<SegyFile*>(this);

  self->hdr_ptr_ =
      first_hdr_ptr_ +
      offset * (sizeof(Trace::Header) + sizeof(float) * num_samples_per_trc_);
  self->trc_ptr_ = hdr_ptr_ + sizeof(Trace::Header);
  self->cur_offset_ = offset;
}

bool SegyFile::read(Trace& trace) const {
  checkFileOpened();

  Trace::Header& header = trace.header();
  std::vector<float>& samples = trace.samples();

  if (hdr_ptr_ >= file_->end() ||
      (hdr_ptr_ + sizeof(Trace::Header)) >= file_->end()) {
    return false;
  }

  std::memcpy(&header, hdr_ptr_, sizeof(Trace::Header));

  size_t trace_bytes = sizeof(float) * num_samples_per_trc_;
  if (trc_ptr_ >= file_->end() || (trc_ptr_ + trace_bytes) > file_->end())
    return false;

  samples.resize(num_samples_per_trc_);
  std::memcpy(&(samples[0]), trc_ptr_, trace_bytes);

  uint16_t sample_format_code =
      binary_header_->getValueAtOffset<uint16_t>(kSegyDataSampleFormat);

  if (sample_format_code == kSegyIEEEDataFormatCode) {
#if __BYTE_ORDER == __LITTLE_ENDIAN
    std::transform(samples.begin(), samples.end(), samples.begin(),
                   [](float a) -> float { return swap_endianness(a); });
#endif
  } else if (sample_format_code == kSegyIBMFloatDataFormatCode) {
    std::transform(samples.begin(), samples.end(), samples.begin(),
                   [](float a) -> float { return ibm_to_ieee(a, true); });
  } else {
    throw std::runtime_error("Segy: read: Data format not supported : " +
                             std::to_string(sample_format_code));
  }

  return true;
}

const SegyFile::TextHeader& SegyFile::getTextHeader() const {
  checkFileOpened();
  return *text_header_;
}

void SegyFile::setTextHeader(const SegyFile::TextHeader& header) {
  checkFileOpened();
  (*text_header_) = header;
}

const SegyFile::BinaryHeader& SegyFile::getBinaryHeader() const {
  checkFileOpened();
  return *binary_header_;
}

void SegyFile::setBinaryHeader(const SegyFile::BinaryHeader& header) {
  checkFileOpened();
  (*binary_header_) = header;
}

std::string SegyFile::TextHeader::toString() const {
  std::ostringstream ostr;
  bool isEBCDIC = isTextInEBCDICFormat(&data[0], sizeof(data));
  std::string line(kSegyTextHeaderLineWidth, ' ');
  for (size_t i = 0, j = 0; i < sizeof(data); i++, j++) {
    j = j % kSegyTextHeaderLineWidth;
    if (isEBCDIC) {
      line[j] = getASCIIForEBCDIC(data[i]);
    } else {
      line[j] = data[i];
    }

    if ((i + 1) % kSegyTextHeaderLineWidth == 0) {
      ostr << line << std::endl;
    }
  }

  return ostr.str();
}

void SegyFile::BinaryHeader::print(std::ostream& os) const {
  os << "SEGY Version : "
     << int(getValueAtOffset<uint8_t>(kSegyBinHdrSegyMajorOffset)) << "."
     << int(getValueAtOffset<uint8_t>(kSegyBinHdrSegyMinorOffset)) << std::endl;
  os << "Sample interval : "
     << getValueAtOffset<uint16_t>(kSegyBinHdrSamplingIntOffset) << std::endl;
  os << "Num samples per trc : "
     << getValueAtOffset<uint16_t>(kSegyNumSamplesPerTrcOffset) << std::endl;
  uint16_t format_code = getValueAtOffset<uint16_t>(kSegyDataSampleFormat);
  os << "Data sample format code : " << format_code << " (";
  switch (format_code) {
    case 1:
      os << "IBM floating-point";
      break;
    case 5:
      os << "IEEE floating-point";
      break;
    default:
      os << "Not supported";
      break;
  }
  os << ")" << std::endl;
  uint16_t units = getValueAtOffset<uint16_t>(55);
  os << "Measurement units: " << ((units == 1) ? "Meters" : "Feet")
     << std::endl;
}

void SegyFile::Trace::Header::print(std::ostream& os) const {
  os << std::endl;

  os << "Possible inline/crossline locations: " << std::endl;
  for (size_t i = 0; i < kSegyHeaderInlineCandidateOffsets.size(); i++) {
    int il_offset = kSegyHeaderInlineCandidateOffsets[i];
    int xl_offset = kSegyHeaderCrosslineCandidateOffsets[i];

    int32_t il = getValueAtOffset<int32_t>(il_offset);
    int32_t xl = getValueAtOffset<int32_t>(xl_offset);
    os << "Offset (" << il_offset << ", " << xl_offset << ") -> (" << il << ", "
       << xl << ")" << std::endl;
  }

  os << std::endl << "Possible coordinate locations: " << std::endl;
  for (size_t i = 0; i < kSegyHeaderXCoordCandidateOffsets.size(); i++) {
    int x_offset = kSegyHeaderXCoordCandidateOffsets[i];
    int y_offset = kSegyHeaderYCoordCandidateOffsets[i];

    float x_coord = getCoordinateValue(x_offset);
    float y_coord = getCoordinateValue(y_offset);
    os << "Offset (" << x_offset << ", " << y_offset << ") -> (" << x_coord
       << ", " << y_coord << ")" << std::endl;
  }
}

void SegyFile::checkFileOpened() const {
  if (!file_->is_open()) {
    std::ostringstream ostr;
    ostr << "File '" << file_->name() << "' has not been opened!" << std::endl;
    throw std::runtime_error(ostr.str());
  }
}

void SegyFile::checkFileOpenedForReading() const {
  if (!file_->is_open() || mode_ & std::ios_base::in) {
    std::ostringstream ostr;
    ostr << "File '" << file_->name() << "' has not been opened!" << std::endl;
    throw std::runtime_error(ostr.str());
  }
}

}  // namespace segystack
