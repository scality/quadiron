/* -*- mode: c++ -*- */
#pragma once

template <typename CharT>
class ostreambuf : public std::basic_streambuf<CharT>
{
public:
  ostreambuf(CharT* buffer, size_t buf_len)
  {
    std::basic_streambuf<CharT>::setp(buffer, buffer + buf_len);
  }
};

template <typename CharT>
class istreambuf : public std::basic_streambuf<CharT>
{
public:
  istreambuf(CharT* buffer, size_t buf_len)
  {
    std::basic_streambuf<CharT>::setg(buffer, buffer, buffer + buf_len);
  }
};
