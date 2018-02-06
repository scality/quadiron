#!/usr/bin/env ruby
# encoding: utf-8

COMPIL_OPTS = '-x c++ --std=c++11 -Isrc'.freeze
CHECKS = %w[
  android-*
  boost-*
  clang-analyzer-*
  mpi-*
].freeze

# These files don't compile for now
BLACKLIST = %w[
  tools/find_pminus1_compo.cpp
  tools/find_primes.cpp
  tools/is_prime.cpp
].freeze
TO_CHECK = Dir.glob('{benchmark,src,test,tools}/*.{cpp,h}') - BLACKLIST

output = `clang-tidy #{TO_CHECK.join(' ')} -checks='#{CHECKS.join(',')}' -- #{COMPIL_OPTS}`
if output =~ /(error|warning):/
  puts output
  exit(1)
end
