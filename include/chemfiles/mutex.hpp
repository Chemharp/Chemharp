// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_MUTEX_HPP
#define CHEMFILES_MUTEX_HPP

#include <mutex>
#include <utility>

namespace chemfiles {

template<class T> class mutex;

/// A lock guard that guarantee exclusive access to the underlying data.
///
/// When the lock guard is destroyed, it releases the associated mutex.
template<class T>
class lock_guard {
public:
    T& operator*() {
        return data_;
    }

    T* operator->() {
        return &data_;
    }

    const T& operator*() const {
        return data_;
    }

    const T* operator->() const {
        return &data_;
    }

private:
    lock_guard(T& data, std::unique_lock<std::mutex>&& guard):
        data_(data), guard_(std::move(guard)) {}

    T& data_;
    std::unique_lock<std::mutex> guard_;

    friend class mutex<T>;
};

/// A `std::mutex` wrapper that own some associated data, and ensure that the
/// mutex is always locked before accessing the data
template<class T>
class mutex {
public:
    /// Create a new mutex containing the given `data`
    mutex(T data): data_(std::move(data)) {}
    /// Create a new mutex containing a default constructed T
    mutex(): mutex(T()) {}

    mutex(const mutex&) = delete;
    mutex& operator=(const mutex&) = delete;
    mutex(mutex&&) = delete;
    mutex& operator=(mutex&&) = delete;

    ~mutex() {
        // deadlock if someone is trying to destroy this `chemfiles::mutex`
        // while still holding a lock on the underlying `mutex_`.
        mutex_.lock();
        mutex_.unlock();
    }

    /// Lock the mutex, and return a `lock_guard`. The `lock_guard` allow to
    /// access the locked data, and will release the mutex when it goes out of
    /// scope.
    lock_guard<T> lock() {
        return lock_guard<T>(data_, std::unique_lock<std::mutex>(mutex_));
    }

private:
    T data_;
    std::mutex mutex_;
};

} // namespace chemfiles

#endif
