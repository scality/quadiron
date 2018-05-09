/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 Scality
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __QUAD_PROPERTY_H__
#define __QUAD_PROPERTY_H__

#include <cstdint>
#include <iosfwd>
#include <string>
#include <unordered_map>

#include <sys/types.h>

namespace quad {

/** The location of a value in a buffer. */
struct ValueLocation {
    off_t offset;
    uint32_t fragment_id;

    /** Construct an invalid location. */
    ValueLocation() : offset(-1), fragment_id(0) {}

    /** Construct a location from an offset and a fragment ID. */
    ValueLocation(off_t pos, uint32_t fragment)
        : offset(pos), fragment_id(fragment)
    {
    }

    /** Construct a location from its string representation. */
    explicit ValueLocation(const std::string& str);

    /**  Checks if the two locations are equal. */
    bool operator==(const ValueLocation& other) const
    {
        return offset == other.offset && fragment_id == other.fragment_id;
    }

    /** Return the string representation of the location. */
    std::string to_string() const
    {
        std::string str(std::to_string(offset));
        str += ':';
        str += std::to_string(fragment_id);
        return str;
    }

    inline off_t get_offset() const
    {
        return offset;
    }
};

} // namespace quad

namespace std {

template <>
struct hash<quad::ValueLocation> {
    std::size_t operator()(const quad::ValueLocation& k) const
    {
        return std::hash<off_t>()(k.offset)
               + std::hash<uint32_t>()(k.fragment_id);
    }
};

} // namespace std

namespace quad {

/** Ancillary data attached to values.
 *
 * A property carries extra-information (whose interpretation is left to the
 * reader) related to a specific value (identified by its location).
 */
class Properties {
  public:
    inline void add(const ValueLocation& loc, const std::string& data)
    {
        props[loc] = data;
    }

    inline const std::string* const get(const ValueLocation& loc) const
    {
        auto it = props.find(loc);
        return it != props.end() ? &it->second : nullptr;
    }

    inline void clear()
    {
        props.clear();
    }

    std::unordered_map<ValueLocation, std::string> const get_map() const
    {
        return props;
    }

  private:
    std::unordered_map<ValueLocation, std::string> props;

    friend std::istream& operator>>(std::istream& is, Properties& props);
    friend std::ostream& operator<<(std::ostream& os, const Properties& props);
};

} // namespace quad

#endif
