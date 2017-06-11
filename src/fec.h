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
    TYPE_1, //Systematic code: take n_data input, generate n_parities outputs
    TYPE_2  //Non-systematic code: take n_data input, generate n_data+n_parities outputs
  };

  FECType fec_type;
  u_int word_size;
  u_int n_data;
  u_int n_parities;

protected:
  GF<T> *gf;

 public:
  FEC(GF<T> *gf, FECType fec_type, u_int word_size, u_int n_data, u_int n_parities);
  virtual int get_n_fragments_required() = 0;
  virtual int get_n_inputs() = 0;
  virtual int get_n_outputs() = 0;
  virtual void encode(std::vector<KeyValue*> props, off_t offset, Vec<T> *output, Vec<T> *words) = 0;
  virtual void decode_add_data(int fragment_index, int row) = 0;
  virtual void decode_add_parities(int fragment_index, int row) = 0;
  virtual void decode_build(void) = 0;
  virtual void decode(std::vector<KeyValue*> props, off_t offset, Vec<T> *output, Vec<T> *words) = 0;

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
 * @param gf 
 * @param word_size in bytes
 * @param n_data 
 * @param n_parities 
 */
template <typename T>
FEC<T>::FEC(GF<T> *gf, FECType fec_type, u_int word_size, u_int n_data, u_int n_parities)
{
  assert(fec_type == TYPE_1 || fec_type == TYPE_2);
  
  this->gf = gf;
  this->fec_type = fec_type;
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
  } else {
    assert(false && "no such size");
  }
  return false;
}

/** 
 * Encode buffers
 * 
 * @param input_data_bufs must be exactly n_data
 * @param output_parities_bufs must be exactly n_parities (set nullptr when not missing/wanted)
 * @param output_parities_props specific properties that the called is supposed to store along with parities
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
    for (int i = 0;i < n_data;i++) {
      T tmp;
      if (!readw(&tmp, input_data_bufs[i])) {
        cont = false;
        break ;
      }
      words.set(i, tmp);
    }
    if (!cont)
      break ;
    encode(output_parities_props, offset, &output, &words);
    for (int i = 0;i < get_n_outputs();i++) {
      T tmp = output.get(i);
      writew(tmp, output_parities_bufs[i]);
    }
    offset += word_size;
  } 
}

/** 
 * Decode buffers
 * 
 * @param input_data_bufs if TYPE_1 must be exactly n_data otherwise 0 (use nullptr when missing)
 * @param input_parities_bufs if TYPE_1 must be exactly n_parities otherwise n_data+n_parities (use nullptr when missing)
 * @param input_parities_props caller is supposed to provide specific information bound to parities
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

  assert(input_data_bufs.size() == n_data);
  assert(input_parities_bufs.size() == get_n_outputs());
  assert(input_parities_props.size() == get_n_outputs());
  assert(output_data_bufs.size() == n_data);

  if (fec_type == TYPE_1) {
    for (int i = 0;i < n_data;i++) {
      if (input_data_bufs[i] != nullptr) {
        decode_add_data(fragment_index, i);
        fragment_index++;
      }
    }

    //data is in clear so nothing to do
    if (fragment_index == n_data)
      return true;
  }
  
  if (fragment_index < get_n_fragments_required()) {
    //finish with parities available
    for (int i = 0;i < get_n_outputs();i++) {
      if (input_parities_bufs[i] != nullptr) {
        decode_add_parities(fragment_index, i);
        fragment_index++;
        //stop when we have enough parities
        if (fragment_index == get_n_fragments_required())
          break ;
      }
    }
    
    //unable to decode
    if (fragment_index < get_n_fragments_required())
      return false;
  }
    
  decode_build();
    
  Vec<T> words(gf, get_n_inputs());
  Vec<T> output(gf, get_n_inputs());
    
  while (true) {
    words.zero_fill();
    fragment_index = 0;
    if (fec_type == TYPE_1) {
      for (int i = 0;i < get_n_inputs();i++) {
        if (input_data_bufs[i] != nullptr) {
          T tmp;
          if (!readw(&tmp, input_data_bufs[i])) {
            cont = false;
            break ;
          }
          words.set(fragment_index, tmp);
          fragment_index++;
        }
      }
      if (!cont)
        break ;
    }
    //stop when we have enough parities
    if (fragment_index == get_n_fragments_required())
      break ;
    for (int i = 0;i < get_n_outputs();i++) {
      if (input_parities_bufs[i] != nullptr) {
        T tmp;
        if (!readw(&tmp, input_parities_bufs[i])) {
          cont = false;
          break ;
        }
        words.set(fragment_index, tmp);
        fragment_index++;
        //stop when we have enough parities
        if (fragment_index == get_n_fragments_required())
          break ;
      }
    }
    if (!cont)
      break ;

    decode(input_parities_props, offset, &output, &words);
      
    for (int i = 0;i < n_data;i++) {
      if (output_data_bufs[i] != nullptr) {
        T tmp = output.get(i);
        writew(tmp, output_data_bufs[i]);
      }
    }

    offset += word_size;
  } 

  return true;
}
