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
#ifndef __QUAD_BENCH_STATS_H__
#define __QUAD_BENCH_STATS_H__

#include <iostream>
#include <string>
#include <cmath>

class Stats_t {
  public:
    Stats_t(const std::string& str, size_t work_load)
    {
        this->name = str;
        this->work_load = work_load;

        this->nb = 0;
        this->sum = 0;
        this->sum_2 = 0;
    }

    void begin()
    {
        nb = 0;
        sum = 0;
        sum_2 = 0;
    }

    void add(uint64_t val)
    {
        nb++;
        sum += val;
        sum_2 += (val * val);
    }

    void end()
    {
        avg = static_cast<double>(sum) / static_cast<double>(nb);
        std_dev = sqrt(
            static_cast<double>(sum_2) / static_cast<double>(nb) - avg * avg);
    }

    void show()
    {
        std::cout << name << ":\tLatency(us) " << avg << " +/- " << std_dev;
        std::cout << "\t\tThroughput " << work_load / avg << " (MB/s)"
                  << std::endl;
    }

    double get_avg()
    {
        return avg;
    }

    double get_std_dev()
    {
        return std_dev;
    }

    double get_thrpt()
    {
        return work_load / avg;
    }

  private:
    uint64_t nb;
    uint64_t sum;
    uint64_t sum_2;
    double avg;
    double std_dev;
    size_t work_load;
    std::string name;
};

#endif
