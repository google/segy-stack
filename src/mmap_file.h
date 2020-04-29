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

#ifndef SEGYSTACK_MMAPFILE_H_
#define SEGYSTACK_MMAPFILE_H_

#include <fstream>
#include <memory>
#include <string>

namespace segystack {

class MmapFile {
 public:
  typedef std::ios_base::openmode OpenMode;

  static std::unique_ptr<MmapFile> Create(const std::string& filename) {
    std::unique_ptr<MmapFile> ptr(new MmapFile(filename));
    return ptr;
  }

  MmapFile(const std::string& file);

  bool exists() const;

  // deletes the file if it exists. file will be closed first if open.
  void remove();

  void open(std::ios_base::openmode mode = std::ios_base::in);

  void close();

  bool is_open() const { return isopen_; }

  const std::string& name() const { return filename_; }

  // Returns the file descriptor associated with this file.
  int fd() const;

  uint64_t size() const { return filesize_; }

  ~MmapFile();

  // Expands the current file by |length| bytes.
  void expand(uint64_t length);

  void flush();

  OpenMode open_mode() const { return mode_; }

  void map();

  void unmap();

  bool is_mapped() const { return is_mapped_; }

  void* addr() { return addr_; }

  char* char_addr() { return reinterpret_cast<char*>(addr_); }

  char* begin() { return char_addr(); }

  char* end() { return (char_addr() + filesize_); }

  std::ptrdiff_t get_addr_offset(char* ptr) const {
    return ptr - reinterpret_cast<char*>(addr_);
  }

 protected:
  MmapFile(const MmapFile&) = delete;
  MmapFile& operator=(const MmapFile&) = delete;

  std::string filename_;
  bool isopen_;
  OpenMode mode_;
  uint64_t filesize_;
  int fd_;
  bool is_mapped_;
  void* addr_;
};

}  // namespace segystack

#endif  // SEGYSTACK_MMAPFILE_H_
