/* -*- mode: c++ -*- */
#ifndef __NTTEC_CONFIG_H__
#define __NTTEC_CONFIG_H__

#include <iosfwd>
#include <map>
#include <string>

struct KeyValue : std::map<std::string, std::string> {
    bool is_key(const std::string& s) const
    {
        return count(s) != 0;
    }
};

std::istream& operator>>(std::istream& ins, KeyValue& d);
std::ostream& operator<<(std::ostream& outs, const KeyValue& d);

#endif
