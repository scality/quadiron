#include <iostream>
#include <sstream>

#include "property.h"

namespace nttec {

ValueLocation::ValueLocation(const std::string& str)
{
    std::istringstream iss(str);
    std::string buf;

    std::getline(iss, buf, ':');
    offset = std::stol(buf);

    std::getline(iss, buf);
    fragment_id = std::stoul(buf);
}

std::istream& operator>>(std::istream& is, Properties& props)
{
    std::string line;

    while (std::getline(is, line)) {
        // Skip blank lines.
        auto begin = line.find_first_not_of(" \f\t\v");
        if (begin == std::string::npos) {
            continue;
        }
        // Skip commentary.
        if (std::string("#;").find(line[begin]) != std::string::npos) {
            continue;
        }
        // Extract the key.
        auto end = line.find('=', begin);
        std::string key = line.substr(begin, end - begin);
        // No leading or trailing whitespace allowed.
        key.erase(key.find_last_not_of(" \f\t\v") + 1);
        // No blank keys allowed
        if (key.empty()) {
            continue;
        }
        // Extract the data (no leading or trailing whitespace allowed).
        begin = line.find_first_not_of(" \f\n\r\t\v", end + 1);
        end = line.find_last_not_of(" \f\n\r\t\v") + 1;

        props.add(ValueLocation(key), line.substr(begin, end - begin));
    }

    return is;
}

std::ostream& operator<<(std::ostream& os, const Properties& props)
{
    for (auto& kv : props.props) {
        os << kv.first.to_string() << " = " << kv.second << '\n';
    }
    return os;
}

} // namespace nttec
