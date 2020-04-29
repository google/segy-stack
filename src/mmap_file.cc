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

#include "mmap_file.h"

#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace segystack {

MmapFile::MmapFile(const std::string& file)
    : filename_(file),
      isopen_(false),
      filesize_(0),
      fd_(-1),
      is_mapped_(false),
      addr_(0) {}

MmapFile::~MmapFile() {
  if (is_mapped_)
    unmap();
  if (isopen_)
    close();
}

bool MmapFile::exists() const {
  struct stat buffer;
  if (::stat(filename_.c_str(), &buffer) != 0) {
    return false;
  }
  return true;
}

void MmapFile::remove() {
  if (is_mapped())
    unmap();
  if (is_open())
    close();
  if (!exists())
    return;

  if (::unlink(filename_.c_str()) != 0) {
    std::cerr << "MmapFile::delete failed on file " << name() << "! "
              << std::endl
              << "Reason: " << ::strerror(errno) << std::endl;
  }
}

void MmapFile::open(std::ios_base::openmode mode) {
  if (isopen_)
    return;

  if (mode & std::ios_base::in) {
    struct stat buffer;
    if (::stat(filename_.c_str(), &buffer) != 0) {
      throw std::runtime_error("File " + filename_ + " does not exist!");
    }
  }

  mode_ = mode;

  std::ostringstream flags_str;
  int open_flags;
  mode_t mode_flags = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH;
  if (mode_ & std::ios_base::in) {
    open_flags = O_RDONLY;
    flags_str << "RDONLY";
  } else {
    open_flags = O_RDWR;
    flags_str << "RDWR";
  }

  if (mode & std::ios_base::trunc) {
    open_flags |= O_CREAT;
    flags_str << "| CREAT";
  }

  if ((fd_ = ::open(filename_.c_str(), open_flags, mode_flags)) < 0) {
    std::ostringstream ostr;
    ostr << "MmapFile::open() - Unable to open file : " << filename_
         << std::endl;
    ostr << "Flags : " << flags_str.str() << std::endl;
    throw std::runtime_error(ostr.str());
  }

  if (mode & std::ios_base::trunc) {
    if (::ftruncate(fd_, 0) != 0) {
      std::cerr << "MmapFile::ftruncate failed on file " << name() << "! "
                << std::endl
                << "Reason: " << ::strerror(errno) << std::endl;
    }
  }

  ::lseek(fd_, 0, SEEK_SET);
  filesize_ = ::lseek(fd_, 0, SEEK_END);

  isopen_ = true;
}

void MmapFile::close() {
  if (!isopen_)
    return;
  if (is_mapped_)
    unmap();

  ::close(fd_);

  isopen_ = false;
}

int MmapFile::fd() const {
  if (isopen_)
    return fd_;
  throw std::runtime_error("File not opened!");
}

void MmapFile::expand(uint64_t length) {
  if (!isopen_ || mode_ & std::ios_base::in) {
    std::ostringstream ostr;
    ostr << "MmapFile(" << filename_
         << "): Cannot expand! File not opened or opened in read only mode!"
         << std::endl;
    throw std::runtime_error(ostr.str());
  }

  uint64_t old_filesize = filesize_;

  int retval = ::ftruncate(fd_, filesize_ + length);

  if (retval != 0) {
    std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << ::strerror(errno)
              << std::endl;
    return;
  }

  filesize_ += length;
  ::lseek(fd_, old_filesize,
          SEEK_SET);  // get back to where we were before expanding the file
}

void MmapFile::flush() {
  if (::msync(addr_, filesize_, MS_SYNC) != 0) {
    std::cerr << "MmapFile::msync failed on file " << name() << "! "
              << std::endl
              << "Reason: " << ::strerror(errno) << std::endl;
  }
  ::fsync(fd_);
}

void MmapFile::map() {
  if (is_mapped_)
    return;

  if (!is_open()) {
    std::ostringstream ostr;
    ostr << "MmapFile(" << name() << "): Cannot map! File not opened!"
         << std::endl;
    throw std::runtime_error(ostr.str());
  }

  int prot_flags, mmap_flags;
  switch (mode_) {
    case std::ios_base::in:
      prot_flags = PROT_READ;
      mmap_flags = MAP_PRIVATE;
      break;
    default:
      prot_flags = PROT_READ | PROT_WRITE;
      mmap_flags = MAP_SHARED;
      break;
  }

  ::lseek(fd_, 0, SEEK_SET);
  filesize_ = ::lseek(fd_, 0, SEEK_END);

  if ((addr_ = ::mmap(NULL, filesize_, prot_flags, mmap_flags, fd_, 0)) ==
      MAP_FAILED) {
    std::ostringstream ostr;
    ostr << "MmapFile::map() - Unable to mmap data for " << name() << std::endl;
    ostr << "Reason: " << ::strerror(errno) << std::endl;
    ostr << "filesize = " << filesize_ << std::endl;
    ostr << "fd = " << fd_ << std::endl;
    throw std::runtime_error(ostr.str());
  }

  is_mapped_ = true;
}

void MmapFile::unmap() {
  if (!is_mapped_)
    return;

  flush();

  if (::munmap(addr_, filesize_) != 0) {
    std::cerr << "MmapFile::munmap failed on file " << name() << "! "
              << std::endl
              << "Reason: " << ::strerror(errno) << std::endl;
  }

  is_mapped_ = false;
}

}  // namespace segystack
