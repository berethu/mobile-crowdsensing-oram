
t_bits my_randseed;

benesFunction t_bits random_bits();


//////
// Auxiliary stuff

typedef enum {right,left} t_direction;

benesFunction t_bool odd(t_int x);
benesFunction t_longint gcd(t_longint a, t_longint b);
benesFunction t_bits mul_inv(t_bits x);
benesFunction t_bits rol(t_bits x, t_int rot);
benesFunction t_bits rol_lo(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits rolc_lo(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits gray_code(t_bits x);
benesFunction t_bits inv_gray_code(t_bits x);
benesFunction t_int nr_1bits(t_bits x);
benesFunction t_int nr_leading_0bits(t_bits x);
benesFunction t_int nr_trailing_0bits(t_bits x);
benesFunction t_bool is_contiguous_1bits(t_bits x);
benesFunction t_bits tbm(t_bits x, t_int mode);

benesFunction t_bits blend(t_bits m, t_bits x, t_bits y);
benesFunction t_bits simd_odd(t_bits x, t_subword sw);

benesFunction t_bits bit_permute_step(t_bits x, t_bits m, t_uint shift);
benesFunction t_bits bit_permute_step_simple(t_bits x, t_bits m, t_uint shift);
benesFunction void bit_permute_step2(t_bits* x1, t_bits* x2, t_bits m, t_uint shift);

benesFunction void identity_perm(ta_index tgt);
benesFunction void invert_perm(const ta_index src, ta_index tgt);
benesFunction void random_perm(ta_index tgt);
benesFunction t_bits used_source_bits(const ta_index perm);
benesFunction t_bits used_target_bits(const ta_index perm);


//////
// Bit index operations

benesFunction t_bits bit_index_complement(t_bits x, t_subword k);
benesFunction t_bits bit_index_swap(t_bits x, t_subword j, t_subword k);
benesFunction t_bits bit_index_swap_complement(t_bits x, t_subword j, t_subword k);

benesFunction t_bits bit_index_ror(t_bits x, t_subword ofs, t_subword field, t_int rot);
benesFunction t_bits transpose(t_bits x, t_subword ld_fields, t_subword ld_col, t_subword ld_row);
benesFunction t_bits shuffle_power(t_bits x, t_subword sw1, t_subword sw2, t_int pwr);
benesFunction t_bits unshuffle_power(t_bits x, t_subword sw1, t_subword sw2, t_int pwr);

benesFunction t_bits permute_bpc(t_bits x, const ta_subword tgt, t_subword_set k);
benesFunction void invert_bpc(const ta_subword src, t_subword_set src_k, ta_subword tgt, t_subword_set* tgt_k);


//////
// Generalized Bit Reversal

benesFunction t_bits general_reverse_bits(t_bits x, t_int k);
benesFunction t_bits bswap(t_bits x);


//////
// Swap by primitives

benesFunction t_bits prim_swap(t_bits x, t_bits m);


//////
// Shuffle and unshuffle

benesFunction t_bits shuffle(t_bits x, t_subword sw1, t_subword sw2);
benesFunction t_bits unshuffle(t_bits x, t_subword sw1, t_subword sw2);


//////
// A "class" for butterfly and other operations

typedef struct {
  // This structure is used to hold the configuration of
  // butterfly-based operations as well as compress and expand.

  t_bits cfg[ld_bits];  // butterfly configuration
  t_bits mask;  // saved mask, for compress/expand

  // Here is sketched how to convert this to a class:
  // Include all the generator and usage benesFunctions as private methods
  // and replace the parameter self by the implicit object pointer this.
  // Add the many compound routines.
  // Remove the name suffix  for all methods.
  // If you want to cache the configuration, add here:
  //   kind: the generator kind
  //     enum (initialized, frot, vrot, ce_right, ce_left, cef_right, cef_left)
  //   sw: the used subword size (t_subword)
  // Add an initializer/constructor which sets kind to initialized.
  // The generator routines must set the keys (kind, mask, sw).
  // The compound routines check the cached keys (kind, mask, sw);
  //   if not equal, call the generator routine and update the configuration;
  // finally they call the usage routine.
  } tr_bfly;


//////
// Compress and expand

//////
// Compress and expand: Compress bit masks

benesFunction t_bits compress_mask_right(t_bits m, t_subword sw);
benesFunction t_bits compress_mask_left(t_bits m, t_subword sw);
benesFunction t_bits compress_mask(t_bits m, t_subword sw, t_direction d);


//////
// Compress and expand: Generate configuration

benesFunction void gen_ce_right(tr_bfly* self, t_bits m, t_subword sw);
benesFunction void gen_ce_left(tr_bfly* self, t_bits m, t_subword sw);


//////
// Compress and expand: Usage

benesFunction t_bits apply_compress_right(const tr_bfly* self, t_bits x);
benesFunction t_bits apply_compress_left(const tr_bfly* self, t_bits x);
benesFunction t_bits apply_expand_right(const tr_bfly* self, t_bits x);
benesFunction t_bits apply_expand_left(const tr_bfly* self, t_bits x);


//////
// Compress and expand: Compound

benesFunction t_bits compress_right(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits compress_left(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits compress(t_bits x, t_bits m, t_subword sw, t_direction d);
benesFunction t_bits expand_right(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits expand_left(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits expand(t_bits x, t_bits m, t_subword sw, t_direction d);


//////
// Butterfly network

benesFunction t_bits butterfly(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits bfly(const tr_bfly* self, t_bits x);
benesFunction t_bits ibfly(const tr_bfly* self, t_bits x);
benesFunction t_bool bfly_parity(const tr_bfly* self);


//////
// Rotate via butterfly

//////
// Rotate via butterfly: Generate configuration

benesFunction void gen_frot(tr_bfly* self, t_int rot, t_subword sw);
benesFunction void gen_vrot(tr_bfly* self, t_bits rot, t_subword sw);


//////
// Rotate via butterfly: Compound

benesFunction t_bits fror_bfly(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits frol_bfly(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits frot_bfly(t_bits x, t_int rot, t_subword sw, t_direction d);

benesFunction t_bits vror_bfly(t_bits x, t_bits rot, t_subword sw);
benesFunction t_bits vrol_bfly(t_bits x, t_bits rot, t_subword sw);
benesFunction t_bits vrot_bfly(t_bits x, t_bits rot, t_subword sw, t_direction d);

benesFunction t_bits frol(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits fror(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits frot(t_bits x, t_int rot, t_subword sw, t_direction d);

benesFunction t_bits frolc(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits frorc(t_bits x, t_int rot, t_subword sw);
benesFunction t_bits frotc(t_bits x, t_int rot, t_subword sw, t_direction d);

benesFunction t_bits vrol(t_bits x, t_bits rot, t_subword sw);
benesFunction t_bits vror(t_bits x, t_bits rot, t_subword sw);
benesFunction t_bits vrot(t_bits x, t_bits rot, t_subword sw, t_direction d);


//////
// Compress/expand-flip via butterfly

//////
// Compress/expand-flip via butterfly: Generate configuration

benesFunction void gen_cef_right(tr_bfly* self, t_bits m, t_subword sw);
benesFunction void gen_cef_left(tr_bfly* self, t_bits m, t_subword sw);


//////
// Compress/expand-flip via butterfly: Compound

benesFunction t_bits compress_flip_right(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits compress_flip_left(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits compress_flip(t_bits x, t_bits m, t_subword sw, t_direction d);
benesFunction t_bits expand_flip_right(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits expand_flip_left(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits expand_flip(t_bits x, t_bits m, t_subword sw, t_direction d);


//////
// Omega/flip

benesFunction t_bits omega(t_bits x, t_bits m, t_subword sw);
benesFunction t_bits flip(t_bits x, t_bits m, t_subword sw);


//////
// Permutations via Benes network

typedef struct {
  tr_bfly b1,b2;
  } tr_benes;

benesFunction void gen_benes_ex(tr_benes* self, const ta_index c_tgt, const ta_subword a_stage);
benesFunction void gen_benes(tr_benes* self, const ta_index c_tgt);
benesFunction t_bits benes_fwd(const tr_benes* self, t_bits x);
benesFunction t_bits benes_bwd(const tr_benes* self, t_bits x);
benesFunction t_bits benes_fwd_ex(const tr_benes* self, t_bits x, const ta_subword a_stage);
benesFunction t_bits benes_bwd_ex(const tr_benes* self, t_bits x, const ta_subword a_stage);
benesFunction t_bool benes_parity(const tr_benes* self);

// eof.
