/* -*- mode: c++ -*- */
#ifndef __NTL_NETL_BENCH_IOSTREAMBUF_H__
#define __NTL_NETL_BENCH_IOSTREAMBUF_H__

#include <ios>
#include <streambuf>

template <typename CharT>
class ostreambuf : public std::basic_streambuf<CharT>
{
 private:
  CharT* begin;
  CharT* end;

 public:
  ostreambuf(CharT* buffer, size_t buf_len)
  {
    this->begin = buffer;
    this->end = buffer + buf_len;
    std::basic_streambuf<CharT>::setp(begin, end);
  }

  std::ios::pos_type seekpos(std::ios::pos_type sp, std::ios_base::openmode which) {
    std::basic_streambuf<CharT>::setp(begin + sp, end);
    std::ios::pos_type rc = std::basic_streambuf<CharT>::seekpos(sp, which);
    return rc;
  }
};

template <typename CharT>
class istreambuf : public std::basic_streambuf<CharT>
{
private:
  CharT* begin;
  CharT* end;

public:
  istreambuf(CharT* buffer, size_t buf_len)
  {
    this->begin = buffer;
    this->end = buffer + buf_len;
    std::basic_streambuf<CharT>::setg(begin, begin, end);
  }

  std::ios::pos_type seekpos(std::ios::pos_type sp, std::ios_base::openmode which) {
    std::basic_streambuf<CharT>::setg(begin + sp, begin + sp, end);
    std::ios::pos_type rc = std::basic_streambuf<CharT>::seekpos(sp, which);
    return rc;
  }
};

#endif
