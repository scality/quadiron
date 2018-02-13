/* -*- mode: c++ -*- */
#ifndef __NTTEC_PROPERTY_H__
#define __NTTEC_PROPERTY_H__

#include <cstdint>
#include <iosfwd>
#include <string>
#include <unordered_map>

#include <sys/types.h>

namespace nttec {

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
};

} // namespace nttec

namespace std {

template <>
struct hash<nttec::ValueLocation> {
    std::size_t operator()(const nttec::ValueLocation& k) const
    {
        return std::hash<off_t>()(k.offset)
               + std::hash<uint32_t>()(k.fragment_id);
    }
};

} // namespace std

namespace nttec {

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

  private:
    std::unordered_map<ValueLocation, std::string> props;

    friend std::istream& operator>>(std::istream& is, Properties& props);
    friend std::ostream& operator<<(std::ostream& os, const Properties& props);
};

} // namespace nttec

#endif
