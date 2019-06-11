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

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iosfwd>
#include <string>
#include <vector>
#include <assert.h>

#include <netinet/in.h>
#include <sys/types.h>

#include "core.h"

namespace quadiron {

static constexpr unsigned OOR_MARK = 1;

/** Ancillary data attached to values.
 *
 * A property carries extra-information (whose interpretation is left to the
 * reader) related to a specific value (identified by its location).
 * It wraps a map whose each element is a key/value where
 *  - key indicates the location of symbol whose value should be adjusted
 *  - value indicates value that could be used to adjust the symbol value
 * For prime fields, value is always 1.
 * For NF4, value is an uint32_t integer.
 */
class Properties {
  private:
    std::vector<std::pair<size_t, uint32_t>> props;
    size_t m_default_key = 0;
    uint32_t m_default_val = 0;

    friend std::istream& operator>>(std::istream& is, Properties& props);
    friend std::ostream& operator<<(std::ostream& os, const Properties& props);

  public:
    enum { FNT1 = 0x464E5431 };

    /**
     * Add a pair of location and marker
     */
    inline void add(const size_t location, const uint32_t marker)
    {
        props.push_back(std::make_pair(location, marker));
    }

    inline void clear()
    {
        props.clear();
    }

    const std::vector<std::pair<size_t, uint32_t>>& get_map() const
    {
        return props;
    }

    /**
     * Sort `props` on the basis of location in ascending order
     */
    inline void sort()
    {
        std::sort(props.begin(), props.end());
    }

    /**
     * Serialize properties into a buffer (FNT)
     *
     * @return 0 if OK, else -1
     */
    inline int fnt_serialize(uint32_t* dwords, unsigned n_dwords)
    {
        if ((2 + props.size()) > n_dwords) {
            return -1;
        }
        dwords[0] = htonl(FNT1);
        unsigned i = 2;
        for (auto const& item : props) {
            dwords[i++] = htonl(narrow_cast<uint32_t>(item.first));
            dwords[i++] = htonl(narrow_cast<uint32_t>(item.second));
        }
        dwords[1] = htonl(i - 2);
        std::fill(dwords + i, dwords + n_dwords - 1, htonl(0));

        return 0;
    }

    /**
     * Deserialize properties from a buffer (FNT)
     *
     * @return 0 if OK, else -1
     */
    inline int fnt_deserialize(const uint32_t* dwords, unsigned n_dwords)
    {
        if (n_dwords < 2) {
            return -1;
        }
        uint32_t magic = ntohl(dwords[0]);
        if (magic != FNT1) {
            return -1;
        }
        uint32_t _n_dwords = ntohl(dwords[1]);
        if ((2 + _n_dwords) > n_dwords) {
            return -1;
        }
        for (unsigned i = 0; i < _n_dwords; i += 2) {
            add(static_cast<size_t>(ntohl(dwords[i + 2])),
                ntohl(dwords[i + 3]));
        }
        return 0;
    }

    /**
     * Get location at a given index
     *
     * @param index - index of element
     *
     * @returns location of the element if the element is found, otherwise a
     * default value is return
     */
    inline const size_t& location(size_t index) const
    {
        return (index < props.size()) ? props[index].first : m_default_key;
    }

    /**
     * Get marker at a given index
     *
     * @param index - index of element
     *
     * @returns marker of the element if the element is found, otherwise a
     * default value is return
     */
    inline const uint32_t& marker(size_t index) const
    {
        return (index < props.size()) ? props[index].second : m_default_val;
    }

    /**
     * Check if a given location is the location of the `index`th element
     *
     * @param index - index of element
     * @param loc - a location
     *
     * @returns true if it matches, false otherwise
     */
    inline bool is_marked(size_t index, const size_t loc) const
    {
        return (index < props.size()) ? (props[index].first == loc) : false;
    }

    /**
     * Check if the location of the `index`th element is within a range
     *
     * @param index - index of element
     * @param min - lower bound of the range (inclusive)
     * @param max - upper bound of the range (exclusive)
     *
     * @returns true if it matches, false otherwise
     */
    inline bool in_range(size_t index, const size_t min, const size_t max) const
    {
        return (index < props.size())
                   ? (props[index].first >= min && props[index].first < max)
                   : false;
    }
};

class FntProperties : public Properties {
  public:
};

} // namespace quadiron

#endif
