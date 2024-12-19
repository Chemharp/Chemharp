// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include "chemfiles/capi/shared_allocator.hpp"
#include "chemfiles/mutex.hpp"

// define the global allocator instance
chemfiles::mutex<chemfiles::shared_allocator> chemfiles::shared_allocator::instance_;
