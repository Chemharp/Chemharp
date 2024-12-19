#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include <memory>
#include <string>
#include <type_traits>

#include "chemfiles/File.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/warnings.hpp"

#include "chemfiles/files/BinaryFile.hpp"

#include <sys/fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <utility>

#ifdef CHEMFILES_WINDOWS
#include <io.h>
#endif

#ifdef __CYGWIN__
    #define fseek64 fseek
    #define ftell64 ftell
    #define off64_t off_t
#elif defined(_MSC_VER)
    #define fseek64 _fseeki64
    #define ftell64 _ftelli64
    #define open _open
    #define fdopen _fdopen
    #define off64_t __int64
#elif !CHEMFILES_BINARY_FILE_USE_MMAP
    // assume unix by default
    #define fseek64 fseeko
    #define ftell64 ftello
    #define off64_t off_t
    static_assert(_FILE_OFFSET_BITS == 64, "_FILE_OFFSET_BITS must be 64");
#endif

#if CHEMFILES_BINARY_FILE_USE_MMAP
    #include <sys/mman.h>
    #include <unistd.h>
#endif

using namespace chemfiles;

BinaryFile::BinaryFile(std::string path, File::Mode mode):
    File(std::move(path), mode, File::Compression::DEFAULT)
{
    int open_mode;
    if (mode == Mode::READ) {
        open_mode = O_RDONLY;
    } else if (mode == Mode::APPEND) {
        open_mode = O_RDWR | O_CREAT;
    } else if (mode == Mode::WRITE) {
        open_mode = O_RDWR | O_CREAT | O_TRUNC;
    } else {
        unreachable();
    }

#ifdef CHEMFILES_WINDOWS
    open_mode |= _O_BINARY;
    int permissions = _S_IWRITE;
#else
    int permissions = S_IRWXU | S_IRWXG | S_IROTH;
#endif

    auto file_descriptor = open(this->path().c_str(), open_mode, permissions);
    if (file_descriptor == -1) {
        throw file_error(
            "could not open file at '{}': {}", this->path(), std::strerror(errno)
        );
    }

#if CHEMFILES_BINARY_FILE_USE_MMAP
    file_descriptor_ = file_descriptor;

    if (mode == Mode::READ) {
        mmap_prot_ = PROT_READ;
    } else {
        mmap_prot_ = PROT_READ | PROT_WRITE;
    }

    if (mode != File::WRITE) {
        struct stat file_stat;
        auto status = fstat(file_descriptor_, &file_stat);
        if (status < 0) {
            throw file_error("could not get the file size with fstat: {}", std::strerror(errno));
        }
        file_size_ = static_cast<size_t>(file_stat.st_size);
    }

    page_size_ = static_cast<size_t>(sysconf(_SC_PAGESIZE));

    if (mode == File::READ) {
        mmap_size_ = file_size_;
        total_written_size_ = file_size_;
    } else if (mode == File::WRITE) {
        // set the map size to a large value (1GiB here) to reduce the number of
        // time we have to unmap/mmap the file. This drastically improves
        // performances when writing, at the cost of only being able to use mmap
        // on 64-bit systems (on 32-bits systems, there is a very low chance to
        // be able to get 1GiB contiguous allocation even in virtual memory, so
        // we fallback to the implementation based on cstdio instead).
        //
        // This does not mean chemfiles will actually hold that much RAM, only a
        // virtual allocation of that size.
        mmap_size_ = 0x40000000;
    } else if (mode == File::APPEND) {
        mmap_size_ = 0x40000000;
        // make sure the file fits inside mmap_size_
        while (file_size_ > mmap_size_) {
            mmap_size_ *= 2;
        }

        offset_ = file_size_;
        total_written_size_ = file_size_;
    } else {
        unreachable();
    }

    mmap_data_ = static_cast<char*>(mmap(
        nullptr, mmap_size_, mmap_prot_, MAP_SHARED, file_descriptor_, 0
    ));

    if (mmap_data_ == MAP_FAILED) {
        throw file_error("mmap failed for '{}': {}", this->path(), std::strerror(errno));
    }
#else
    const char* fdopen_mode;
    if (mode == Mode::READ) {
        fdopen_mode = "rb";
    } else {
        fdopen_mode = "r+b";
    }

    file_ = fdopen(file_descriptor, fdopen_mode);
    assert(file_ != nullptr);

    if (mode == Mode::APPEND) {
        auto status = fseek64(file_, 0, SEEK_END);
        if (status != 0) {
            throw file_error("failed to seek in file: {}", std::strerror(errno));
        }
    } else {
        this->seek(0);
    }
#endif
}

BinaryFile::~BinaryFile() noexcept {
    this->close_file();
}

BinaryFile& BinaryFile::operator=(BinaryFile&& other) noexcept {
    File::operator=(std::move(other));

    this->close_file();

#if CHEMFILES_BINARY_FILE_USE_MMAP
    std::swap(this->file_descriptor_, other.file_descriptor_);
    std::swap(this->total_written_size_, other.total_written_size_);
    std::swap(this->mmap_data_, other.mmap_data_);
    std::swap(this->file_size_, other.file_size_);
    std::swap(this->mmap_size_, other.mmap_size_);
    std::swap(this->mmap_prot_, other.mmap_prot_);
    std::swap(this->page_size_, other.page_size_);
    std::swap(this->offset_, other.offset_);
#else
    std::swap(this->file_, other.file_);
#endif

    return *this;
}

void BinaryFile::close_file() noexcept {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    if (mmap_data_ != nullptr) {
        auto status = msync(mmap_data_, mmap_size_, MS_SYNC);
        if (status != 0) {
            warning(
                "binary file writer",
                "failed to sync file ({}), some data might be lost",
                std::strerror(errno)
            );
        }

        status = munmap(mmap_data_, mmap_size_);
        if (status != 0) {
            warning(
                "binary file writer",
                "failed to unmap file ({}), something might be wrong",
                std::strerror(errno)
            );
        }
    }

    if (file_descriptor_ != -1) {
        if (this->mode() != Mode::READ) {
            auto status = ftruncate(file_descriptor_, static_cast<off_t>(total_written_size_));
            if (status != 0) {
                warning(
                    "binary file writer",
                    "failed to resize file to it's actual size when closing: {}",
                    std::strerror(errno)
                );
            }
        }

        if (close(file_descriptor_) != 0) {
            warning(
                "binary file writer",
                "failed to close the file ({}), something might be wrong",
                std::strerror(errno)
            );
        }
    }

    file_descriptor_ = -1;
    total_written_size_ = 0;
    mmap_data_ = nullptr;
    file_size_ = 0;
    mmap_size_ = 0;
    mmap_prot_ = 0;
    page_size_ = 0;
    offset_ = 0;
#else
    if (file_ == nullptr) {
        return;
    }

    if (this->mode() != Mode::READ) {
        if (std::fflush(file_) != 0) {
            warning("binary file writer", "failed to flush when closing the file, some data might be lost");
        }
    }

    if (std::fclose(file_) != 0) {
        warning("binary file writer", "failed to close the file, something might be wrong");
    }
    file_ = nullptr;
#endif
}

void BinaryFile::read_char(char* data, size_t count) {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    if (offset_ + count > file_size_) {
        throw file_error(
            "failed to read {} bytes from the file at '{}': mmap out of bounds",
            count, this->path()
        );
    }
    std::memcpy(data, mmap_data_ + offset_, count);
    offset_ += count;
#else
    auto read = std::fread(data, 1, count, file_);
    const char* error_info = "unknown cause";
    if (read != count) {
        if (feof(file_)) {
            error_info = "reached end of file";
        }
        else if (ferror(file_)) {
            error_info = std::strerror(errno);
        }
        throw file_error(
            "failed to read {} bytes from the file at '{}': {}",
            count, this->path(), error_info
        );
    }
#endif
}


void BinaryFile::write_char(const char* data, size_t count) {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    if (offset_ + count > file_size_) {
        while (offset_ + count > file_size_) {
            // increase the file size by multiples of the page size to have to
            // call ftruncate less often
            file_size_ += 4 * page_size_;
        }
        // resize the file, but keep the same mmap allocation unless
        // file_size_ > mmap_size_
        auto status = ftruncate(file_descriptor_, static_cast<off_t>(file_size_));
        if (status != 0) {
            throw file_error("failed to resize file: {}", std::strerror(errno));
        }

        if (file_size_ > mmap_size_) {
            // unmap & remap file with bigger mapping
            status = msync(mmap_data_, mmap_size_, MS_SYNC);
            if (status != 0) {
                throw file_error(
                    "failed to sync file ({}), some data might be lost",
                    std::strerror(errno)
                );
            }

            status = munmap(mmap_data_, mmap_size_);
            if (status != 0) {
                throw file_error("failed to unmap file: {}", std::strerror(errno));
            }

            while (file_size_ > mmap_size_) {
                mmap_size_ *= 2;
            }

            mmap_data_ = static_cast<char*>(mmap(
                nullptr, mmap_size_, mmap_prot_, MAP_SHARED, file_descriptor_, 0
            ));

            if (mmap_data_ == MAP_FAILED) {
                throw file_error("mmap failed for '{}': {}", this->path(), std::strerror(errno));
            }
        }
    }

    if (offset_ + count > total_written_size_) {
        total_written_size_ = offset_ + count;
    }

    std::memcpy(mmap_data_ + offset_, data, count);
    offset_ += count;
#else
    auto written = std::fwrite(data, 1, count, file_);
    if (written != count) {
        throw file_error(
            "failed to write {} bytes to the file at '{}': {}",
            count, this->path(), std::strerror(errno)
        );
    }
#endif
}


uint64_t BinaryFile::tell() const {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    return offset_;
#else
    auto position = ftell64(file_);
    if (position < 0) {
        throw file_error("call to ftell failed: {}", std::strerror(errno));
    }
    return static_cast<uint64_t>(position);
#endif
}


void BinaryFile::seek(uint64_t position) {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    offset_ = position;
#else
    auto status = fseek64(file_, static_cast<off64_t>(position), SEEK_SET);

    if (status != 0) {
        throw file_error("failed to seek in file: {}", std::strerror(errno));
    }
#endif
}


void BinaryFile::skip(uint64_t count) {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    offset_ += count;
#else
    auto status = fseek64(file_, static_cast<off64_t>(count), SEEK_CUR);
    if (status != 0) {
        throw file_error("failed to seek in file: {}", std::strerror(errno));
    }
#endif
}


uint64_t BinaryFile::file_size() {
#if CHEMFILES_BINARY_FILE_USE_MMAP
    return total_written_size_;
#else
    auto current = this->tell();

    auto status = fseek64(file_, 0L, SEEK_END);
    if (status != 0) {
        throw file_error("failed to seek in file: {}", std::strerror(errno));
    }

    auto file_size = this->tell();
    this->seek(current);
    return file_size;
#endif
}


/******************************************************************************/

#define CHEMFILES_LITTLE_ENDIAN 0
#define CHEMFILES_BIG_ENDIAN 1

/***** adapted from boost ******/
// GNU libc offers the helpful header <endian.h> which defines
// __BYTE_ORDER
#if defined (__GLIBC__)
# include <endian.h>
# if (__BYTE_ORDER == __LITTLE_ENDIAN)
#  define CHEMFILES_BYTE_ORDER CHEMFILES_LITTLE_ENDIAN
# elif (__BYTE_ORDER == __BIG_ENDIAN)
#  define CHEMFILES_BYTE_ORDER CHEMFILES_BIG_ENDIAN
# else
#  error Unknown machine endianness detected.
# endif
#elif (defined(_BIG_ENDIAN) || defined(__BIG_ENDIAN__)) && !(defined(_LITTLE_ENDIAN) || defined(__LITTLE_ENDIAN__))
# define CHEMFILES_BYTE_ORDER CHEMFILES_BIG_ENDIAN
#elif (defined(_LITTLE_ENDIAN) || defined(__LITTLE_ENDIAN__)) && !(defined(_BIG_ENDIAN) || defined(__BIG_ENDIAN__))
# define CHEMFILES_BYTE_ORDER CHEMFILES_LITTLE_ENDIAN
#elif defined(__sparc) || defined(__sparc__) \
   || defined(_POWER) || defined(__powerpc__) \
   || defined(__ppc__) || defined(__hpux) || defined(__hppa) \
   || defined(_MIPSEB) || defined(_POWER) \
   || defined(__s390__)
# define CHEMFILES_BYTE_ORDER CHEMFILES_BIG_ENDIAN
#elif defined(__i386__) || defined(__alpha__) \
   || defined(__ia64) || defined(__ia64__) \
   || defined(_M_IX86) || defined(_M_IA64) \
   || defined(_M_ALPHA) || defined(__amd64) \
   || defined(__amd64__) || defined(_M_AMD64) \
   || defined(__x86_64) || defined(__x86_64__) \
   || defined(_M_X64) || defined(__bfin__)
# define CHEMFILES_BYTE_ORDER CHEMFILES_LITTLE_ENDIAN
#elif defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
# define CHEMFILES_BYTE_ORDER CHEMFILES_LITTLE_ENDIAN
#elif defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
# define CHEMFILES_BYTE_ORDER CHEMFILES_BIG_ENDIAN
#else
# error Unknown target machine endianness
#endif
/***** adapted from boost ******/

template<size_t Size>
struct unsigned_type;

template<>
struct unsigned_type<2>{
    typedef uint16_t type;
};

template<>
struct unsigned_type<4>{
    typedef uint32_t type;
};

template<>
struct unsigned_type<8>{
    typedef uint64_t type;
};

template<typename T>
inline T swap_endianness(T value) {
    // swap endianness of a value by byte-copying it into an unsigned integer
    // with same size swapping the integer and byte-copying it back.
    // The compilers should see right through this and optimize away the memcpy
    static_assert(std::is_arithmetic<T>::value, "type must be arithmetic for swapping");
    typename unsigned_type<sizeof(T)>::type uval;
    std::memcpy(&uval, &value, sizeof(T));
    uval = swap_endianness(uval);
    std::memcpy(&value, &uval, sizeof(T));
    return value;
}

template<>
inline uint16_t swap_endianness(uint16_t value) {
    return static_cast<uint16_t>(
        ((static_cast<uint32_t>(value) & 0xff) << 8) |
        ((static_cast<uint32_t>(value) >> 8) & 0xff)
    );
}

template<>
inline uint32_t swap_endianness(uint32_t value) {
    return (
        (value << 24) |
        ((value <<  8) & 0x00ff0000) |
        ((value >>  8) & 0x0000ff00) |
        ((value >> 24) & 0x000000ff)
    );
}

template<>
inline uint64_t swap_endianness(uint64_t value) {
    return (
        ((value & 0x00000000000000FFULL) << 56) |
        ((value & 0x000000000000FF00ULL) << 40) |
        ((value & 0x0000000000FF0000ULL) << 24) |
        ((value & 0x00000000FF000000ULL) <<  8) |
        ((value & 0x000000FF00000000ULL) >>  8) |
        ((value & 0x0000FF0000000000ULL) >> 24) |
        ((value & 0x00FF000000000000ULL) >> 40) |
        ((value & 0xFF00000000000000ULL) >> 56)
    );
}

/******************************************************************************/

template<typename T>
inline void BigEndianFile::read_as_big_endian(T* data, size_t count) {
    auto* char_data = reinterpret_cast<char*>(data);
    const size_t byte_count = sizeof(T) * count;
    this->read_char(char_data, byte_count);
#if CHEMFILES_BYTE_ORDER == CHEMFILES_LITTLE_ENDIAN
    for (size_t i = 0; i < count; ++i) {
        char* value_ptr = char_data + i * sizeof(T);
        T value;
        std::memcpy(&value, value_ptr, sizeof(T));
        value = swap_endianness(value);
        std::memcpy(value_ptr, &value, sizeof(T));
    }
#endif
}

template<typename T>
inline void BigEndianFile::write_as_big_endian(const T* data, size_t count) {
    const size_t byte_count = sizeof(T) * count;
#if CHEMFILES_BYTE_ORDER == CHEMFILES_LITTLE_ENDIAN
    swap_buffer_.resize(byte_count);
    for (size_t i = 0; i < count; ++i) {
        T value = swap_endianness(data[i]);
        std::memcpy(swap_buffer_.data() + i * sizeof(T), &value, sizeof(T));
    }
    this->write_char(swap_buffer_.data(), byte_count);
#else
    this->write_char(reinterpret_cast<const char*>(data), byte_count);
#endif
}

void BigEndianFile::read_i16(int16_t* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_u16(uint16_t* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_i32(int32_t* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_u32(uint32_t* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_i64(int64_t* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_u64(uint64_t* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_f32(float* data, size_t count) {
    read_as_big_endian(data, count);
}


void BigEndianFile::read_f64(double* data, size_t count) {
    read_as_big_endian(data, count);
}

void BigEndianFile::write_i16(const int16_t* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_u16(const uint16_t* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_i32(const int32_t* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_u32(const uint32_t* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_i64(const int64_t* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_u64(const uint64_t* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_f32(const float* data, size_t count) {
    write_as_big_endian(data, count);
}

void BigEndianFile::write_f64(const double* data, size_t count) {
    write_as_big_endian(data, count);
}


/******************************************************************************/

template<typename T>
inline void LittleEndianFile::read_as_little_endian(T* data, size_t count) {
    auto* char_data = reinterpret_cast<char*>(data);
    const size_t byte_count = sizeof(T) * count;
    this->read_char(char_data, byte_count);
#if CHEMFILES_BYTE_ORDER == CHEMFILES_BIG_ENDIAN
    for (size_t i = 0; i < count; ++i) {
        char* value_ptr = char_data + i * sizeof(T);
        T value;
        std::memcpy(&value, value_ptr, sizeof(T));
        value = swap_endianness(value);
        std::memcpy(value_ptr, &value, sizeof(T));
    }
#endif
}

template<typename T>
inline void LittleEndianFile::write_as_little_endian(const T* data, size_t count) {
    const size_t byte_count = sizeof(T) * count;
#if CHEMFILES_BYTE_ORDER == CHEMFILES_BIG_ENDIAN
    swap_buffer_.resize(byte_count);
    for (size_t i = 0; i < count; ++i) {
        T value = swap_endianness(data[i]);
        std::memcpy(swap_buffer_.data() + i * sizeof(T), &value, sizeof(T));
    }
    this->write_char(swap_buffer_.data(), byte_count);
#else
    this->write_char(reinterpret_cast<const char*>(data), byte_count);
#endif
}


void LittleEndianFile::read_i16(int16_t* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_u16(uint16_t* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_i32(int32_t* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_u32(uint32_t* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_i64(int64_t* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_u64(uint64_t* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_f32(float* data, size_t count) {
    read_as_little_endian(data, count);
}


void LittleEndianFile::read_f64(double* data, size_t count) {
    read_as_little_endian(data, count);
}

void LittleEndianFile::write_i16(const int16_t* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_u16(const uint16_t* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_i32(const int32_t* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_u32(const uint32_t* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_i64(const int64_t* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_u64(const uint64_t* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_f32(const float* data, size_t count) {
    write_as_little_endian(data, count);
}

void LittleEndianFile::write_f64(const double* data, size_t count) {
    write_as_little_endian(data, count);
}

/******************************************************************************/

std::unique_ptr<BinaryFile> BinaryFile::open_native(std::string path, File::Mode mode) {
#if CHEMFILES_BYTE_ORDER == CHEMFILES_LITTLE_ENDIAN
    return std::make_unique<LittleEndianFile>(std::move(path), mode);
#else
    return std::make_unique<BigEndianFile>(std::move(path), mode);
#endif
}
