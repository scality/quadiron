/* -*- mode: c++ -*- */
#ifndef __NTTEC_FEC_RS_GFP_FFT_H__
#define __NTTEC_FEC_RS_GFP_FFT_H__

#include "arith.h"
#include "fec_base.h"
#include "fft_2n.h"
#include "fft_base.h"
#include "fft_ct.h"
#include "gf_prime.h"
#include "polynomial.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace nttec {
namespace fec {

/** Reed-Solomon (RS) Erasure code over prime Galois Fields and FFT.
 *
 * @note Necessary condition: p ≥ 2<sup>8*word_size</sup> to cover all
 * possible values of symbols, as each symbol is 8*word_size bits.
 *
 * A length-k vector X is encoded into a length-n vector Y whose elements are
 * GF(p) elements. Hence, elements of Y could be greater than
 * 2<sup>8*word_size</sup>
 *
 * We can deal with this as below.
 *
 * @note To facilitate our implementation, we choose
 *   p < 2 * 2<sup>8*word_size</sup>.
 *
 * Supposing that an element Y[i] > 2<sup>8*word_size</sup>, we will store its
 * value modulo 2<sup>8*word_size</sup>, i.e.
 *
 *   Y<sub>p</sub>[i] = Y[i] % 2<sup>8*word_size</sup>
 *
 * A flag will be used to mark the position `i`.
 *
 * To decode X, we read Y[i] from Y<sub>p</sub>[i] as:
 *
 *    Y[i] = Y<sub>p</sub>[i] + 2<sup>8*word_size</sup>
 *
 * Because p < 2 * 2<sup>8*word_size</sup>, a single bool is enough as flag.
 */
template <typename T>
class RsGfpFft : public FecCode<T> {
  public:
    T n;
    T r;

    RsGfpFft(unsigned word_size, unsigned n_data, unsigned n_parities)
        : FecCode<T>(FecType::NON_SYSTEMATIC, word_size, n_data, n_parities)
    {
        // warning all fermat numbers >= to F_5 (2^32+1) are composite!!!
        T gf_p = 0;
        if (word_size < 4) {
            gf_p = (1ULL << (8 * word_size)) + 1;
            this->limit_value = (1ULL << (8 * word_size));
        } else if (word_size == 4) {
            gf_p = (T)4294991873ULL;              // p-1=2^13 29^1 101^1 179^1
            this->limit_value = (T)4294967296ULL; // 2^32
        } else {
            assert(false); // not support yet
            exit(1);
        }

        assert(gf_p >= this->limit_value);
        // we choose gf_p for a simple implementation
        assert(gf_p / 2 < this->limit_value);

        this->gf = new gf::Prime<T>(gf_p);
        assert(
            arith::jacobi<T>(this->gf->get_primitive_root(), this->gf->card())
            == -1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = this->gf->get_code_len_high_compo(n_parities + n_data);

        // compute root of order n-1 such as r^(n-1) mod q == 1
        r = this->gf->get_nth_root(n);

        // std::cerr << "limit_value=" << limit_value << "\n";
        // std::cerr << "gf_p=" << gf_p << "\n";
        // std::cerr << "n=" << n << "\n";
        // std::cerr << "r=" << r << "\n";

        if (arith::is_power_of_2<T>(n)) {
            this->fft = new fft::Radix2<T>(this->gf, n);
        } else {
            this->fft = new fft::CooleyTukey<T>(this->gf, n);
        }
    }

    ~RsGfpFft()
    {
        delete fft;
        delete this->gf;
    }

    int get_n_outputs()
    {
        return this->n;
    }

    /**
     * Encode vector
     *
     * @param output must be n
     * @param props must be exactly n
     * @param offset used to locate special values
     * @param words must be n_data
     */
    void encode(
        vec::Vector<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words)
    {
        vec::ZeroExtended<T> vwords(words, n);
        fft->fft(output, &vwords);
        // check for out of range value in output
        for (unsigned i = 0; i < this->code_len; i++) {
            if (output->get(i) >= this->limit_value) {
                props[i].add(ValueLocation(offset, i), "@");
                output->set(i, output->get(i) % this->limit_value);
            }
        }
    }

    void decode_add_data(int fragment_index, int row)
    {
        // not applicable
        assert(false);
    }

    void decode_add_parities(int fragment_index, int row)
    {
        // we can't anticipate here
    }

    void decode_build()
    {
        // nothing to do
    }

    /**
     * Perform a Lagrange interpolation to find the coefficients of the
     * polynomial
     *
     * @note If all fragments are available ifft(words) is enough
     *
     * @param output must be exactly n_data
     * @param props special values dictionary must be exactly n_data
     * @param offset used to locate special values
     * @param fragments_ids unused
     * @param words v=(v_0, v_1, ..., v_k-1) k must be exactly n_data
     */
    void decode(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words)
    {
        int k = this->n_data; // number of fragments received
        // vector x=(x_0, x_1, ..., x_k-1)
        vec::Vector<T> vx(this->gf, k);
        for (int i = 0; i < k; i++) {
            vx.set(i, this->gf->exp(r, fragments_ids->get(i)));
        }

        for (int i = 0; i < k; i++) {
            const int j = fragments_ids->get(i);
            auto data = props[j].get(ValueLocation(offset, j));

            if (data && *data == "@") {
                words->set(i, words->get(i) + this->limit_value);
            }
        }

        // Lagrange interpolation
        Polynomial<T> A(this->gf), _A(this->gf);

        // compute A(x) = prod_j(x-x_j)
        A.set(0, 1);
        for (int i = 0; i < k; i++) {
            Polynomial<T> _t(this->gf);
            _t.set(1, 1);
            _t.set(0, this->gf->sub(0, vx.get(i)));
            // _t.dump();
            A.mul(&_t);
        }
        // std::cout << "A(x)="; A.dump();

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        _A.copy(&A);
        _A.derivative();
        // std::cout << "A'(x)="; _A.dump();

        // evaluate n_i=v_i/A'_i(x_i)
        vec::Vector<T> _n(this->gf, k);
        for (int i = 0; i < k; i++) {
            _n.set(i, this->gf->div(words->get(i), _A.eval(vx.get(i))));
        }
        // std::cout << "_n="; _n.dump();

        // compute N'(x) = sum_i{n_i * x^z_i}
        Polynomial<T> N_p(this->gf);
        for (int i = 0; i <= k - 1; i++) {
            N_p.set(fragments_ids->get(i), _n.get(i));
        }

        // We have to find the numerator of the following expression:
        // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
        // using Taylor series we rewrite the expression into
        // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
        Polynomial<T> S(this->gf);
        for (T i = 0; i <= n - 1; i++) {
            T val = this->gf->inv(this->gf->exp(r, i + 1));
            S.set(i, N_p.eval(val));
        }
        S.neg();
        S.mul(&A);
        // No need to mod x^n since only last n_data coefs are obtained
        // output is n_data length
        for (unsigned i = 0; i < this->n_data; i++)
            output->set(i, S.get(i));
    }

  private:
    fft::FourierTransform<T>* fft = nullptr;
    T limit_value;
};

} // namespace fec
} // namespace nttec

#endif
