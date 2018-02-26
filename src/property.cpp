/*
 * Copyright 2017-2018 the NTTEC authors
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
