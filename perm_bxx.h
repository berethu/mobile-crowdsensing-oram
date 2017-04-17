// Local include file to define some derived constants.

#define benesBits (1 << ld_bits)
#define lo_bit (t_bits)(1)
#define hi_bit (lo_bit << (benesBits-1))
#define all_bits (t_bits)(-1)

typedef t_int t_bit_index;  // subrange 0..benesBits-1; -1:don't care
  // used to specify bit indexes

typedef t_bit_index t_subword_set;
  // used to specify a set of t_subword

#define no_index ((t_bit_index)(-1))
  // don't care / wildcard

typedef t_uint t_subword;  // subrange 0..ld_bits
  // used to specify log_2 of subword size, see a_subword

typedef t_subword ta_subword[ld_bits];
  // bit index indexes

typedef t_bit_index ta_index[benesBits];
  // bit indexes

const t_char* a_subword[]={
             // sw benesBits
  "Bit",     // 0  1
  "Nyp",     // 1  2
  "Nibble",  // 2  4
  "Byte",    // 3  8
  "Word",    // 4  16
  "DWord",   // 5  32
  "QWord",   // 6  64
  "OWord",   // 7  128
  "YWord",   // 8  256
  "ZWord"};  // 9  512

