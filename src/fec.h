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
  virtual u_int get_n_outputs() = 0;
  virtual void encode(Vec<T> *output, Vec<T> *words) = 0;
  virtual void repair_init(void) = 0;
  virtual void repair_add_data(int k, int i) = 0;
  virtual void repair_add_parities(int k, int i) = 0;
  virtual void repair_build(void) = 0;
  virtual void repair(Vec<T> *output, Vec<T> *words) = 0;

  bool readw(T *ptr, std::istream *stream);  
  bool writew(T val, std::ostream *stream);  

  void encode_bufs(std::vector<std::istream*> input_data_bufs,
                   std::vector<std::ostream*> output_parities_bufs);
  
  bool decode_bufs(std::vector<std::istream*> input_data_bufs,
                   std::vector<std::istream*> input_parities_bufs,
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
 *
 * @note all streams must be of equal size
 */
template <typename T>
void FEC<T>::encode_bufs(std::vector<std::istream*> input_data_bufs,
                         std::vector<std::ostream*> output_bufs)
{
  bool cont = true;

  assert(input_data_bufs.size() == n_data);
  assert(output_bufs.size() == get_n_outputs());
 
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
    encode(&output, &words);
    for (int i = 0;i < get_n_outputs();i++) {
      T tmp = output.get(i);
      writew(tmp, output_bufs[i]);
    }
  } 
}

/** 
 * Decode channel  
 * 
 * @param input_data_bufs if TYPE_1 must be exactly n_data otherwise 0 (use nullptr when missing)
 * @param input_parities_bufs if TYPE_1 must be exactly n_parities otherwise n_data+n_parities (use nullptr when missing)
 * @param output_data_bufs must be exactly n_data (use nullptr when not missing/wanted)
 *
 * @note All streams must be of equal size
 *
 * @return true if decode succeded, else false
 */
template <typename T>
bool FEC<T>::decode_bufs(std::vector<std::istream*> input_data_bufs,
                         std::vector<std::istream*> input_parities_bufs,
                         std::vector<std::ostream*> output_data_bufs)
{
  bool cont = true;

  u_int n_data_ok = 0;
  u_int n_parities_ok = 0;

  assert(input_data_bufs.size() == n_data);
  assert(input_parities_bufs.size() == get_n_outputs());
  assert(output_data_bufs.size() == n_data);

  repair_init();

  int k = 0;
  if (fec_type == TYPE_1) {
    for (int i = 0;i < n_data;i++) {
      if (input_data_bufs[i] != nullptr) {
        repair_add_data(k, i);
        k++;
        n_data_ok++;
      }
    }
  }

  //XXX optim
  //if (k == n_data)
  //return true;
  
  if (k < n_data) {
    //finish with parities available
    for (int i = 0;i < get_n_outputs();i++) {
      if (input_parities_bufs[i] != nullptr) {
        repair_add_parities(k, i);
        k++;
        n_parities_ok++;
        //stop when we have enough parities
        if (n_data == k)
          break ;
      }
    }
    
    if (n_parities_ok < (n_data - n_data_ok))
      return false;
  }
    
  repair_build();
    
  //read-and-repair
  Vec<T> words(gf, n_data);
  Vec<T> output(gf, n_data);
    
  while (true) {
    words.zero_fill();
    k = 0;
    if (fec_type == TYPE_1) {
      for (int i = 0;i < n_data;i++) {
        if (input_data_bufs[i] != nullptr) {
          T tmp;
          if (!readw(&tmp, input_data_bufs[i])) {
            cont = false;
            break ;
          }
          words.set(k, tmp);
          k++;
        }
      }
      if (!cont)
        break ;
    }
    //stop when we have enough parities
    if (n_data == k)
      break ;
    for (int i = 0;i < get_n_outputs();i++) {
      if (input_parities_bufs[i] != nullptr) {
        T tmp;
        if (!readw(&tmp, input_parities_bufs[i])) {
          cont = false;
          break ;
        }
        words.set(k, tmp);
        k++;
        //stop when we have enough parities
        if (n_data == k)
          break ;
      }
    }
    if (!cont)
      break ;

    repair(&output, &words);
      
    for (int i = 0;i < n_data;i++) {
      if (output_data_bufs[i] != nullptr) {
        T tmp = output.get(i);
        writew(tmp, output_data_bufs[i]);
      }
    }
  } 

  return true;
}
