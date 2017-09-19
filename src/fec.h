/* -*- mode: c++ -*- */
#pragma once

/**
 * Generic class for Forward-Encoding-Codes
 */
template<typename T>
class FEC
{
 public:
  enum FECType
  {
    TYPE_1,  // Systematic code: take n_data input, generate n_parities outputs
    TYPE_2   // Non-systematic code: take n_data input,
             // generate n_data+n_parities outputs
  };

  FECType type;
  u_int word_size;
  u_int n_data;
  u_int n_parities;

  uint64_t total_encode_cycles = 0;
  uint64_t n_encode_ops = 0;
  uint64_t total_decode_cycles = 0;
  uint64_t n_decode_ops = 0;

 protected:
  GF<T> *gf;

 public:
  FEC(FECType type, u_int word_size, u_int n_data, u_int n_parities);

  /**
   * Return the number actual parities for TYPE_1 it is exactly n_parities, for
   * TYPE_2 it maybe at least n_data+n_parities (but sometimes more).
   * @return
   */
  virtual int get_n_outputs() = 0;
  virtual void encode(Vec<T> *output, std::vector<KeyValue*> props,
    off_t offset, Vec<T> *words) = 0;
  virtual void decode_add_data(int fragment_index, int row) = 0;
  virtual void decode_add_parities(int fragment_index, int row) = 0;
  virtual void decode_build(void) = 0;
  /**
   * Decode a vector of words
   *
   * @param props properties bound to parity fragments
   * @param offset offset in the data fragments
   * @param output original data (must be of n_data length)
   * @param fragments_ids identifiers of available fragments
   * @param words input words if TYPE_1 must be n_data,
   *  if TYPE_2 get_n_outputs()
   */
  virtual void decode(Vec<T> *output, std::vector<KeyValue*> props,
    off_t offset, Vec<T> *fragments_ids, Vec<T> *words) = 0;

  bool readw(T *ptr, std::istream *stream);
  bool writew(T val, std::ostream *stream);

  void encode_bufs(std::vector<std::istream*> input_data_bufs,
                   std::vector<std::ostream*> output_parities_bufs,
                   std::vector<KeyValue*> output_parities_props);

  bool decode_bufs(std::vector<std::istream*> input_data_bufs,
                   std::vector<std::istream*> input_parities_bufs,
                   std::vector<KeyValue*> input_parities_props,
                   std::vector<std::ostream*> output_data_bufs);
};

/**
 * Create an encoder
 *
 * @param word_size in bytes
 * @param n_data
 * @param n_parities
 */
template <typename T>
FEC<T>::FEC(FECType type, u_int word_size, u_int n_data, u_int n_parities)
{
  assert(type == TYPE_1 || type == TYPE_2);

  this->type = type;
  this->word_size = word_size;
  this->n_data = n_data;
  this->n_parities = n_parities;
}

template <typename T>
bool FEC<T>::readw(T *ptr, std::istream *stream)
{
  if (word_size == 1) {
    u_char c;
    if (stream->read((char *) &c, 1)) {
      *ptr = c;
      return true;
    }
  } else if (word_size == 2) {
    u_short s;
    if (stream->read((char *) &s, 2)) {
      *ptr = s;
      return true;
    }
  } else if (word_size == 4) {
    u_int s;
    if (stream->read((char *) &s, 4)) {
      *ptr = s;
      return true;
    }
  } else if (word_size == 8) {
    u_long s;
    if (stream->read((char *) &s, 8)) {
      *ptr = s;
      return true;
    }
  } else if (word_size == 16) {
    __uint128_t s;
    if (stream->read((char *) &s, 16)) {
      *ptr = s;
      return true;
    }
  } else {
    assert(false && "no such size");
  }
  return false;
}

template <typename T>
bool FEC<T>::writew(T val, std::ostream *stream)
{
  if (word_size == 1) {
    u_char c = val;
    if (stream->write((char *) &c, 1))
      return true;
  } else if (word_size == 2) {
    u_short s = val;
    if (stream->write((char *) &s, 2))
      return true;
  } else if (word_size == 4) {
    u_int s = val;
    if (stream->write((char *) &s, 4))
      return true;
  } else if (word_size == 8) {
    u_long s = val;
    if (stream->write((char *) &s, 8))
      return true;
  } else if (word_size == 16) {
    __uint128_t s = val;
    if (stream->write((char *) &s, 16))
      return true;
  } else {
    assert(false && "no such size");
  }
  return false;
}

/**
 * Encode buffers
 *
 * @param input_data_bufs must be exactly n_data
 * @param output_parities_bufs must be exactly get_n_outputs() (set nullptr when not missing/wanted)
 * @param output_parities_props must be exactly get_n_outputs() specific properties that the called is supposed to store along with parities
 *
 * @note all streams must be of equal size
 */
template <typename T>
void FEC<T>::encode_bufs(std::vector<std::istream*> input_data_bufs,
                         std::vector<std::ostream*> output_parities_bufs,
                         std::vector<KeyValue*> output_parities_props)
{
  bool cont = true;
  off_t offset = 0;

  assert(input_data_bufs.size() == n_data);
  assert(output_parities_bufs.size() == get_n_outputs());
  assert(output_parities_props.size() == get_n_outputs());

  Vec<T> words = Vec<T>(gf, n_data);
  Vec<T> output = Vec<T>(gf, get_n_outputs());

  while (true) {
    words.zero_fill();
    for (int i = 0; i < n_data; i++) {
      T tmp;
      if (!readw(&tmp, input_data_bufs[i])) {
        cont = false;
        break;
      }
      words.set(i, tmp);
    }
    if (!cont)
      break;

    uint64_t start = rdtsc();
    encode(&output, output_parities_props, offset, &words);
    uint64_t end = rdtsc();
    total_encode_cycles += end - start;
    n_encode_ops++;

    for (int i = 0; i < get_n_outputs(); i++) {
      T tmp = output.get(i);
      writew(tmp, output_parities_bufs[i]);
    }
    offset += word_size;
  }
}

/**
 * Decode buffers
 *
 * @param input_data_bufs if TYPE_1 must be exactly n_data otherwise it is unused (use nullptr when missing)
 * @param input_parities_bufs if TYPE_1 must be exactly n_parities otherwise get_n_outputs() (use nullptr when missing)
 * @param input_parities_props if TYPE_1 must be exactly n_parities otherwise get_n_outputs() caller is supposed to provide specific information bound to parities
 * @param output_data_bufs must be exactly n_data (use nullptr when not missing/wanted)
 *
 * @note All streams must be of equal size
 *
 * @return true if decode succeded, else false
 */
template <typename T>
bool FEC<T>::decode_bufs(std::vector<std::istream*> input_data_bufs,
                         std::vector<std::istream*> input_parities_bufs,
                         std::vector<KeyValue*> input_parities_props,
                         std::vector<std::ostream*> output_data_bufs)
{
  off_t offset = 0;
  bool cont = true;

  u_int fragment_index = 0;

  if (type == TYPE_1)
    assert(input_data_bufs.size() == n_data);
  assert(input_parities_bufs.size() == get_n_outputs());
  assert(input_parities_props.size() == get_n_outputs());
  assert(output_data_bufs.size() == n_data);

  if (type == TYPE_1) {
    for (int i = 0; i < n_data; i++) {
      if (input_data_bufs[i] != nullptr) {
        decode_add_data(fragment_index, i);
        fragment_index++;
      }
    }
    // data is in clear so nothing to do
    if (fragment_index == n_data)
      return true;
  }

  if (fragment_index < n_data) {
    // finish with parities available
    for (int i = 0; i < get_n_outputs(); i++) {
      if (input_parities_bufs[i] != nullptr) {
        decode_add_parities(fragment_index, i);
        fragment_index++;
        // stop when we have enough parities
        if (fragment_index == n_data)
          break;
      }
    }
    // unable to decode
    if (fragment_index < n_data)
      return false;
  }

  decode_build();

  int n_words;
  if (type == TYPE_1)
    n_words = n_data;
  else if (type == TYPE_2)
    n_words = get_n_outputs();

  Vec<T> words(gf, n_words);
  Vec<T> fragments_ids(gf, n_words);
  Vec<T> output(gf, n_data);

  while (true) {
    words.zero_fill();
    fragment_index = 0;
    if (type == TYPE_1) {
      for (int i = 0; i < n_data; i++) {
        if (input_data_bufs[i] != nullptr) {
          T tmp;
          if (!readw(&tmp, input_data_bufs[i])) {
            cont = false;
            break;
          }
          fragments_ids.set(fragment_index, i);
          words.set(fragment_index, tmp);
          fragment_index++;
        }
      }
      if (!cont)
        break;
      // stop when we have enough parities
      if (fragment_index == n_data)
        break;
    }
    for (int i = 0; i < get_n_outputs(); i++) {
      if (input_parities_bufs[i] != nullptr) {
        T tmp;
        if (!readw(&tmp, input_parities_bufs[i])) {
          cont = false;
          break;
        }
        fragments_ids.set(fragment_index, i);
        words.set(fragment_index, tmp);
        fragment_index++;
        // stop when we have enough parities
        if (fragment_index == n_data)
          break;
      }
    }
    if (!cont)
      break;

    uint64_t start = rdtsc();
    decode(&output, input_parities_props, offset, &fragments_ids, &words);
    uint64_t end = rdtsc();
    total_decode_cycles += end - start;
    n_decode_ops++;

    for (int i = 0; i < n_data; i++) {
      if (output_data_bufs[i] != nullptr) {
        T tmp = output.get(i);
        writew(tmp, output_data_bufs[i]);
      }
    }

    offset += word_size;
  }

  return true;
}
