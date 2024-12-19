#ifndef CHFL_SHARED_ALLOCATOR_HPP
#define CHFL_SHARED_ALLOCATOR_HPP

#include <cstddef>
#include <vector>
#include <functional>
#include <unordered_map>
#include <type_traits>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/mutex.hpp"

namespace chemfiles {

struct shared_metadata {
    /// Number of pointer sharing this reference. No need to use atomic
    /// references counting, as the allocator is protected by a mutex
    long count;
    /// How to delete the pointer when we are done with it
    std::function<void(void)> deleter;
};

namespace {
    [[noreturn]] inline void UNINITIALIZED_DELETER() {
        throw chemfiles::error("internal error: uninitialized deleter called");
    }
}

/// An allocator with shared_ptr like semantics, working with raw pointers.
///
/// This is used in the C API to ensure that when taking pointers to
/// atoms/residues/cell inside a frame/topology, the frame/topology is kept
/// alive even if the user calls chfl_free.
class shared_allocator {
public:
    shared_allocator() = default;
    shared_allocator(const shared_allocator&) = delete;
    shared_allocator& operator=(const shared_allocator&) = delete;
    shared_allocator(shared_allocator&&) = default;
    shared_allocator& operator=(shared_allocator&&) = default;

    /// Like `std::make_shared`: create a new shared pointer by constructing a
    /// value of type T with the given arguments.
    template<class T, typename ... Args, typename std::enable_if<!std::is_array<T>::value>::type* = nullptr>
    static T* make_shared(Args&& ... args) {
        auto instance = instance_.lock();
        auto ptr = new T{std::forward<Args>(args)...};
        instance->insert_new(ptr);
        return ptr;
    }

    /// Like `std::make_shared`: create a new shared pointer to an array type.
    /// This function returns a pointer to the first element of the array.
    template<class T, typename std::enable_if<std::is_array<T>::value>::type* = nullptr>
    static typename std::remove_extent<T>::type* make_shared(size_t count) {
        auto instance = instance_.lock();
        auto ptr = new typename std::remove_extent<T>::type[count];
        instance->insert_new_array(ptr);
        return ptr;
    }

    /// Like `std::shared_ptr<U>` aliasing contructor: element and ptr will
    /// share the references count, and none will be freed while the other one
    /// is alive.
    ///
    /// `ptr` must have been allocated with make_shared.
    template<class T, class U>
    static T* shared_ptr(U* ptr, T* element) {
        instance_.lock()->insert_shared(ptr, element);
        return element;
    }

    template<class T, class U>
    static const T* shared_ptr(U* ptr, const T* element) {
        // const_cast is OK here, since we return a const T* anyway, and only
        // use the address of the pointer to create a new shared_ptr
        return shared_ptr(ptr, const_cast<T*>(element));
    }

    /// Decrease the reference count of `ptr`, and delete it if needed.
    static void free(const void* ptr) {
        if (ptr == nullptr) {
            return;
        }
        auto instance = instance_.lock();
        instance->release(ptr);
    }

private:
    template<class T>
    void insert_new(T* ptr) {
        if (pointers_.count(ptr) != 0) {
            throw chemfiles::memory_error(
                "internal error: pointer at {} is already managed by "
                "shared_allocator", static_cast<void*>(ptr)
            );
        }
        size_t id = get_unused_metadata();
        metadata_[id] = shared_metadata{1, [ptr](){ delete ptr; }};
        pointers_.emplace(ptr, id);
    }

    template<class T>
    void insert_new_array(T* ptr) {
        if (pointers_.count(ptr) != 0) {
            throw chemfiles::memory_error(
                "internal error: pointer at {} is already managed by "
                "shared_allocator", static_cast<void*>(ptr)
            );
        }
        size_t id = get_unused_metadata();
        metadata_[id] = shared_metadata{1, [ptr](){ delete[] ptr; }};
        pointers_.emplace(ptr, id);
    }

    void insert_shared(const void* ptr, void* element) {
        auto it = pointers_.find(ptr);
        if (it == pointers_.end()) {
            // the main pointer is not a shared pointer
            throw chemfiles::memory_error(
                "internal error: pointer at {} is not managed by "
                "shared_allocator", ptr
            );
        }

        if (pointers_.count(element) != 0) {
            // Make sure all instances of element in the multi map share the
            // same main pointer & metadata block
            auto id = pointers_.find(element)->second;
            if (id != it->second) {
                // the element pointer is already registered, but with a
                // different main pointer
                throw chemfiles::memory_error(
                    "internal error: element pointer at {} is already managed by "
                    "shared_allocator (associated with {})", element, ptr
                );
            }
        }

        // Insert the new shared pointer
        pointers_.emplace(element, it->second);
        metadata_.at(it->second).count++;
    }

    void release(const void* ptr) {
        auto it = pointers_.find(ptr);
        if (it == pointers_.end()) {
            throw chemfiles::memory_error(
                "unknown pointer passed to shared_allocator::free: {}", ptr
            );
        }

        // Extract the metadat id. We can not just use it->second everywhere,
        // as it can become invalid after the call to `pointers_.erase` below.
        auto id = it->second;
        if (id >= metadata_.size()) {
            throw chemfiles::memory_error(
                "internal error: metadata index is too big: {} >= {}", id, metadata_.size()
            );
        }

        // Decrease refcount
        metadata_[id].count--;
        // Delete the pointer from the pointers map, do not run destructor yet
        pointers_.erase(it);

        if (metadata_[id].count == 0) {
            // Run the destructor and release memory
            metadata_[id].deleter();
            // Mark the metadata block for reuse.
            metadata_[id].deleter = UNINITIALIZED_DELETER;
            unused_.emplace_back(id);
        } else if (metadata_[id].count < 0) {
            throw chemfiles::memory_error(
                "internal error: negative reference count for {}", ptr
            );
        }
    }

    shared_metadata& metadata(const void* ptr) {
        auto it = pointers_.find(ptr);
        if (it != pointers_.end()) {
            return metadata_.at(it->second);
        } else {
            throw chemfiles::memory_error(
                "internal error: unknown pointer passed to shared_allocator::metadata"
            );
        }
    }

    size_t get_unused_metadata() {
        if (!unused_.empty()) {
            // Get an existing one
            size_t id = unused_.back();
            unused_.pop_back();
            return id;
        } else {
            // create a new one
            metadata_.push_back({0, UNINITIALIZED_DELETER});
            return metadata_.size() - 1;
        }
    }

    /// A map of pointer addresses -> indexes of metadata in metadata_
    std::unordered_multimap<const void*, size_t> pointers_;
    /// Metadata for all known pointers
    std::vector<shared_metadata> metadata_;
    /// unused indexes in metadata_ that can be re-used. This is set by
    /// free and used by get_unused_count.
    std::vector<size_t> unused_;

    /// Global instance of the allocator
    static mutex<shared_allocator> instance_;
};

}

#endif
