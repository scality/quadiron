
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "ntl.h"

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class GF2N<uint64_t>;
template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class FEC<uint32_t>;
template class FECGF2NRS<uint32_t>;
template class FECGF2NRS<uint64_t>;
template class FECFNTRS<uint32_t>;
template class FECGF2NFFTRS<uint32_t>;
template class FECGF2NFFTRS<uint64_t>;
template class FECGF2NFFTADDRS<uint32_t>;
template class FECGF2NFFTADDRS<uint64_t>;
template class FECGF2NFFTADDRS<__uint128_t>;
template class FECGFPFFTRS<uint64_t>;

int vflag = 0;
int tflag = 0;
char* prefix = nullptr;

void xusage()
{
    std::cerr << std::string("Usage: ") +
    "ec [-e gf2nrsv|gf2nrsc|gf2nfftrs|gf2nfftaddrs|gfpfftrs|fntrs|ngff4rs]" +
    "[-w word_size][-n n_data][-m n_parities][-p prefix][-v (verbose)]" +
    " -c (encode) | -r (repair)\n";
    exit(1);
}

void xperror(const char* str)
{
    std::cerr << str << "\n";
    exit(1);
}

char* xstrdup(const char* str)
{
    char* n;

    if (nullptr == (n = strdup(str)))
        xperror("malloc");
    return n;
}

/**
 * (re-)create missing prefix.c1 ... cm files
 *
 */
template <typename T>
void create_coding_files(FEC<T>* fec, bool operation_on_packet = false)
{
    char filename[1024];
    std::vector<std::istream*> d_files(fec->n_data, nullptr);
    std::vector<std::ostream*> c_files(fec->n_outputs, nullptr);
    std::vector<std::ostream*> c_props_files(fec->n_outputs, nullptr);
    std::vector<KeyValue*> c_props(fec->n_outputs, nullptr);

    for (unsigned i = 0; i < fec->n_data; i++) {
        snprintf(filename, sizeof(filename), "%s.d%d", prefix, i);
        if (vflag)
            std::cerr << "create: opening data " << filename << "\n";
        d_files[i] = new std::ifstream(filename);
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        snprintf(filename, sizeof(filename), "%s.c%d", prefix, i);
        if (vflag)
            std::cerr << "create: opening coding for writing " << filename
                      << "\n";
        c_files[i] = new std::ofstream(filename);
        snprintf(filename, sizeof(filename), "%s.c%d.props", prefix, i);
        if (vflag)
            std::cerr << "create: opening coding props for writing " << filename
                      << "\n";
        c_props_files[i] = new std::ofstream(filename);
        c_props[i] = new KeyValue();
    }

    if (operation_on_packet)
        fec->encode_packet(d_files, c_files, c_props);
    else
        fec->encode_bufs(d_files, c_files, c_props);

    for (unsigned i = 0; i < fec->n_data; i++) {
        (static_cast<std::ifstream*>(d_files[i]))->close();
        delete d_files[i];
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        *(c_props_files[i]) << *(c_props[i]);

        (static_cast<std::ofstream*>(c_props_files[i]))->close();
        delete c_props_files[i];
        delete c_props[i];

        (static_cast<std::ofstream*>(c_files[i]))->close();
        delete c_files[i];
    }
}

/**
 * repair data files
 *
 */
template <typename T>
bool repair_data_files(FEC<T>* fec)
{
    char filename[1024];
    std::vector<std::istream*> d_files(fec->n_data, nullptr);
    std::vector<std::istream*> c_files(fec->n_outputs, nullptr);
    std::vector<std::istream*> c_props_files(fec->n_outputs, nullptr);
    std::vector<KeyValue*> c_props(fec->n_outputs, nullptr);
    std::vector<std::ostream*> r_files(fec->n_data, nullptr);

    // re-read data
    for (unsigned i = 0; i < fec->n_data; i++) {
        snprintf(filename, sizeof(filename), "%s.d%d", prefix, i);
        if (vflag)
            std::cerr << "repair: checking data " << filename << "\n";
        if (-1 == access(filename, F_OK)) {
            if (vflag)
                std::cerr << filename << " is missing\n";
            d_files[i] = nullptr;
            r_files[i] = new std::ofstream(filename);
        } else {
            r_files[i] = nullptr;
            d_files[i] = new std::ifstream(filename);
        }
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        snprintf(filename, sizeof(filename), "%s.c%d", prefix, i);
        if (vflag)
            std::cerr << "repair: checking coding " << filename << "\n";
        if (access(filename, F_OK)) {
            if (vflag)
                std::cerr << filename << " is missing\n";
            c_files[i] = nullptr;
        } else {
            c_files[i] = new std::ifstream(filename);
        }

        snprintf(filename, sizeof(filename), "%s.c%d.props", prefix, i);
        if (vflag)
            std::cerr << "repair: checking coding props " << filename << "\n";
        if (access(filename, F_OK)) {
            c_props_files[i] = nullptr;
            c_props[i] = nullptr;
        } else {
            c_props_files[i] = new std::ifstream(filename);
            c_props[i] = new KeyValue();
            *(c_props_files[i]) >> *(c_props[i]);
        }
    }

    fec->decode_bufs(d_files, c_files, c_props, r_files);

    for (unsigned i = 0; i < fec->n_data; i++) {
        if (nullptr != d_files[i]) {
            (static_cast<std::ifstream*>(d_files[i]))->close();
            delete d_files[i];
        }
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        if (nullptr != c_props_files[i]) {
            (static_cast<std::ifstream*>(c_props_files[i]))->close();
            delete c_props_files[i];
        }

        if (nullptr != c_props[i])
            delete c_props[i];

        if (nullptr != c_files[i]) {
            (static_cast<std::ifstream*>(c_files[i]))->close();
            delete c_files[i];
        }
    }

    for (unsigned i = 0; i < fec->n_data; i++) {
        if (nullptr != r_files[i]) {
            (static_cast<std::ofstream*>(r_files[i]))->close();
            delete r_files[i];
        }
    }

    return 0;
}

template <typename T>
void print_stats(FEC<T>* fec)
{
    std::cerr << "enc,"
              << (fec->n_encode_ops != 0
                      ? fec->total_encode_cycles / fec->n_encode_ops
                      : 0)
              << ",";
    std::cerr << "dec,"
              << (fec->n_decode_ops != 0
                      ? fec->total_decode_cycles / fec->n_decode_ops
                      : 0)
              << ",";
}

template <typename T>
void print_fec_type(FEC<T>* fec)
{
    switch (fec->type) {
    case FEC<T>::TYPE_1:
        std::cout << "type_1\n";
        break;
    case FEC<T>::TYPE_2:
        std::cout << "type_2\n";
        break;
    default:
        std::cout << "unknown\n";
        break;
    }
}

enum gf2nrs_type {
    VANDERMONDE = 0,
    CAUCHY,
};

template <typename T>
void run_FECGF2NRS(
    int word_size,
    int n_data,
    int n_parities,
    gf2nrs_type mflag,
    int rflag)
{
    FECGF2NRS<T>* fec;
    typename FECGF2NRS<T>::FECGF2NRSType gf2nrs_type;
    ;
    if (mflag == VANDERMONDE)
        gf2nrs_type = FECGF2NRS<T>::VANDERMONDE;
    else
        gf2nrs_type = FECGF2NRS<T>::CAUCHY;
    fec = new FECGF2NRS<T>(word_size, n_data, n_parities, gf2nrs_type);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_FECGF2NFFTRS(int word_size, int n_data, int n_parities, int rflag)
{
    FECGF2NFFTRS<T>* fec;
    fec = new FECGF2NFFTRS<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_FECGF2NFFTADDRS(int word_size, int n_data, int n_parities, int rflag)
{
    FECGF2NFFTADDRS<T>* fec;
    fec = new FECGF2NFFTADDRS<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_FECGFPFFTRS(unsigned word_size, int n_data, int n_parities, int rflag)
{
    assert(sizeof(T) > word_size);

    FECGFPFFTRS<T>* fec;
    fec = new FECGFPFFTRS<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_FECFNTRS(int word_size, int n_data, int n_parities, int rflag)
{
    FECFNTRS<T>* fec;
    size_t pkt_size = 1024;
    fec = new FECFNTRS<T>(word_size, n_data, n_parities, pkt_size);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec, true);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_FECNGFF4RS(int word_size, int n_data, int n_parities, int rflag)
{
    FECNGFF4RS<T>* fec;
    fec = new FECNGFF4RS<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

enum ec_type {
    EC_TYPE_UNDEF = 0,
    EC_TYPE_GF2NRS,
    EC_TYPE_GF2NFFTRS,
    EC_TYPE_GF2NFFTADDRS,
    EC_TYPE_GFPFFTRS,
    EC_TYPE_FNTRS,
    EC_TYPE_NGFF4RS,
};

bool check(int n, int word_size, ec_type eflag)
{
    // we suppose that code length is not too long, i.e. > 2^32
    if (word_size >= 4)
        return true;
    if (eflag == EC_TYPE_FNTRS) {
        return (n <= (1LL << (8 * word_size)) + 1);
    } else {
        return (n <= (1LL << (8 * word_size)));
    }
}

int main(int argc, char** argv)
{
    int n_data, n_parities, opt;
    int cflag = 0;
    int rflag = 0;
    int uflag = 0;
    ec_type eflag = EC_TYPE_UNDEF;
    gf2nrs_type mflag = VANDERMONDE;
    unsigned word_size = 0;

    n_data = n_parities = -1;
    while ((opt = getopt(argc, argv, "n:m:p:cruve:w:t")) != -1) {
        switch (opt) {
        case 'e':
            if (!strcmp(optarg, "gf2nrsv")) {
                eflag = EC_TYPE_GF2NRS;
                mflag = VANDERMONDE;
            } else if (!strcmp(optarg, "gf2nrsc")) {
                eflag = EC_TYPE_GF2NRS;
                mflag = CAUCHY;
            } else if (!strcmp(optarg, "gf2nfftrs")) {
                eflag = EC_TYPE_GF2NFFTRS;
            } else if (!strcmp(optarg, "gf2nfftaddrs")) {
                eflag = EC_TYPE_GF2NFFTADDRS;
            } else if (!strcmp(optarg, "gfpfftrs")) {
                eflag = EC_TYPE_GFPFFTRS;
            } else if (!strcmp(optarg, "ngff4rs")) {
                eflag = EC_TYPE_NGFF4RS;
            } else if (!strcmp(optarg, "fntrs"))
                eflag = EC_TYPE_FNTRS;
            else
                xusage();
            break;
        case 'w':
            word_size = atoi(optarg);
            break;
        case 'v':
            vflag = 1;
            break;
        case 'u':
            uflag = 1;
            break;
        case 'c':
            cflag = 1;
            break;
        case 'r':
            rflag = 1;
            break;
        case 'n':
            n_data = atoi(optarg);
            break;
        case 'm':
            n_parities = atoi(optarg);
            break;
        case 'p':
            prefix = xstrdup(optarg);
            break;
        case 't':
            tflag = 1;
            break;
        default: /* '?' */
            xusage();
        }
    }

    if (!eflag)
        xusage();

    if (!(uflag || cflag || rflag))
        xusage();

    if (0 == check(n_data + n_parities, word_size, eflag)) {
        std::cerr
            << "Number of fragments is too big compared to Galois field size\n";
        exit(1);
    }

    if (-1 == n_data || -1 == n_parities || nullptr == prefix)
        xusage();

    if (eflag == EC_TYPE_FNTRS) {
        if (word_size <= 4)
            run_FECFNTRS<uint32_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 8)
            run_FECFNTRS<uint64_t>(word_size, n_data, n_parities, rflag);
        // else if (word_size <= 16)
        // run_FECFNTRS<__uint128_t>(word_size, n_data, n_parities, rflag);
    } else if (eflag == EC_TYPE_NGFF4RS) {
        if (word_size <= 2)
            run_FECNGFF4RS<uint32_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 4)
            run_FECNGFF4RS<uint64_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 8)
            run_FECNGFF4RS<__uint128_t>(word_size, n_data, n_parities, rflag);
    } else if (eflag == EC_TYPE_GF2NRS) {
        if (word_size <= 4)
            run_FECGF2NRS<uint32_t>(
                word_size, n_data, n_parities, mflag, rflag);
        else if (word_size <= 8)
            run_FECGF2NRS<uint64_t>(
                word_size, n_data, n_parities, mflag, rflag);
        else if (word_size <= 16)
            run_FECGF2NRS<__uint128_t>(
                word_size, n_data, n_parities, mflag, rflag);
    } else if (eflag == EC_TYPE_GFPFFTRS) {
        if (word_size <= 7)
            run_FECGFPFFTRS<uint64_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 15)
            run_FECGFPFFTRS<__uint128_t>(word_size, n_data, n_parities, rflag);
    } else if (eflag == EC_TYPE_GF2NFFTRS) {
        if (word_size <= 4)
            run_FECGF2NFFTRS<uint32_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 8)
            run_FECGF2NFFTRS<uint64_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 16)
            run_FECGF2NFFTRS<__uint128_t>(word_size, n_data, n_parities, rflag);
    } else if (eflag == EC_TYPE_GF2NFFTADDRS) {
        if (word_size <= 4)
            run_FECGF2NFFTADDRS<uint32_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 8)
            run_FECGF2NFFTADDRS<uint64_t>(word_size, n_data, n_parities, rflag);
        else if (word_size <= 16)
            run_FECGF2NFFTADDRS<__uint128_t>(
                word_size, n_data, n_parities, rflag);
    }

    free(prefix);
    return 0;
}
