/* Chemfiles, an efficient IO library for chemistry file formats
* Copyright (C) 2015 Guillaume Fraux
*
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#ifndef CHEMFILES_ERROR_HPP
#define CHEMFILES_ERROR_HPP
#include <stdexcept>
#include <string>

#include "chemfiles/exports.hpp"

namespace chemfiles {

//! @exception Error Base exception class for chemfiles library
struct CHFL_EXPORT Error : public std::runtime_error {
    Error(const std::string& err) : std::runtime_error(err) {}
};

//! @exception FileError Exception for files related failures
struct CHFL_EXPORT FileError : public Error {
    FileError(const std::string& err) : Error(err) {}
};

//! @exception MemoryError Exception for memory related failures
struct CHFL_EXPORT MemoryError : public Error {
    MemoryError(const std::string& err) : Error(err) {}
};

//! @exception FormatError Exception for formats related failures
struct CHFL_EXPORT FormatError : public Error {
    FormatError(const std::string& err) : Error(err) {}
};

//! @exception PluginError Exception for dynamic library loading errors
struct CHFL_EXPORT PluginError : public Error {
    PluginError(const std::string& err) : Error(err) {}
};

//! @exception LexerError Exception for syntaxic errors in selections
class LexerError: public Error {
public:
    LexerError(const std::string& err): Error(err) {}
};

//! @exception ParserError Exception for semantic and parsing errors in selections
class ParserError: public Error {
public:
    ParserError(const std::string& err): Error(err) {}
};

} // namespace chemfiles

#endif
