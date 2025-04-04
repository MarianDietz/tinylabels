#pragma once

#include "seal/seal.h"
#include "seal/util/clipnormal.h"
#include "seal/util/iterator.h"
#include "seal/util/polyarithsmallmod.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;
using namespace seal;
using namespace seal::util;

const size_t mod_plaintext = 50;
const size_t mod_noise = 59;

const size_t poly_modulus_degree = 4096;
const size_t w = 512;  // needs to be a power of 2
const size_t l = 9;    // = log_2 w
const size_t m = 4;    // s.t. g^m > modulus
const uint64_t g = (uint64_t)1 << 28;

const size_t poly_size = 2*poly_modulus_degree;

constexpr double noise_small_standard_deviation = 4;
constexpr double noise_small_max_deviation = 128 * noise_small_standard_deviation;

constexpr double noise_large_standard_deviation = 1000;
constexpr double noise_large_max_deviation = 128 * noise_large_standard_deviation;

void sample_poly_normal(
    shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination,
    double noise_standard_deviation, double noise_max_deviation);

void add_poly_error(
    size_t count,
    shared_ptr<UniformRandomGenerator> prng, const SEALContext::ContextData &context_data, uint64_t *destination,
    double noise_standard_deviation, double noise_max_deviation);

void sample_poly_uniform(
    shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination);

string time_str(chrono::nanoseconds time);
void print_statistics();

struct LHE {
public:

    LHE(const SEALContext::ContextData &context_data, shared_ptr<UniformRandomGenerator> prng) : context_data_(context_data), prng(prng) {}

    void setup();
    void save_pp(FILE* f) {
        fwrite(data_a_.get(), 8, w*poly_size, f);
    }
    void read_pp(FILE* f) {
        data_a_ = allocate_poly_array(w, poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_a_.get(), 8, w*poly_size, f);
    }

    void enc1(Pointer<uint64_t> &m1);
    void save_st1(FILE* f) {
        fwrite(data_s1_.get(), 8, m*poly_size, f);
    }
    void read_st1(FILE* f) {
        data_s1_ = allocate_poly_array(m, poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_s1_.get(), 8, m*poly_size, f);
    }
    void save_ct1(FILE* f) {
        fwrite(data_ct1_.get(), 8, w*m*poly_size, f);
    }
    void read_ct1(FILE* f) {
        data_ct1_ = allocate_poly_array(w*m, poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_ct1_.get(), 8, w*m*poly_size, f);
    }

    void enc2(Pointer<uint64_t> &m2);
    void save_st2(FILE* f) {
        fwrite(data_s2_.get(), 8, poly_size, f);
    }
    void read_st2(FILE* f) {
        data_s2_ = allocate_poly(poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_s2_.get(), 8, poly_size, f);
    }
    void save_ct2(FILE* f) {
        fwrite(data_ct2_.get(), 8, w*poly_size, f);
    }
    void read_ct2(FILE* f) {
        data_ct2_ = allocate_poly_array(w, poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_ct2_.get(), 8, w*poly_size, f);
    }

    void keygen(Pointer<uint64_t> &y);
    void save_sk(FILE* f) {
        fwrite(data_sk_.get(), 8, poly_size, f);
    }
    void read_sk(FILE* f) {
        data_sk_ = allocate_poly(poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_sk_.get(), 8, poly_size, f);
    }

    Pointer<uint64_t>& dec(Pointer<uint64_t> &y);

//private:
    const SEALContext::ContextData &context_data_;
    shared_ptr<UniformRandomGenerator> prng;

    Pointer<uint64_t> data_a_;

    Pointer<uint64_t> data_s1_;
    Pointer<uint64_t> data_s2_;
    Pointer<uint64_t> data_sk_;

    Pointer<uint64_t> data_ct1_;
    Pointer<uint64_t> data_ct2_;

    Pointer<uint64_t> data_mres_;
};

struct Lenc {
public:

    Lenc(const SEALContext::ContextData &context_data, shared_ptr<UniformRandomGenerator> prng) : context_data_(context_data), prng(prng) {}

    void setup();
    void save_pp(FILE* f) {
        fwrite(data_b_.get(), 8, 2*m*poly_size, f);
    }
    void read_pp(FILE* f) {
        data_b_ = allocate_poly_array(2*m, poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_b_.get(), 8, 2*m*poly_size, f);
    }

    Pointer<uint64_t>& enc(Pointer<uint64_t> &s);
    void save_ct1(FILE* f) {
        fwrite(data_ct_.get(), 8, l*w*2*m*poly_size, f);
    }
    void read_ct1(FILE* f) {
        data_ct_ = allocate_poly_array(l*w*2*m, poly_modulus_degree, context_data_.parms().coeff_modulus().size(), MemoryManager::GetPool());
        fread(data_ct_.get(), 8, l*w*2*m*poly_size, f);
    }

    Pointer<uint64_t>& digest(Pointer<uint64_t> &a);

    Pointer<uint64_t>& eval(Pointer<uint64_t> &a);

//private:
    const SEALContext::ContextData &context_data_;
    shared_ptr<UniformRandomGenerator> prng;

    Pointer<uint64_t> data_b_;

    Pointer<uint64_t> data_r_;

    Pointer<uint64_t> data_ct_;

    Pointer<uint64_t> data_tree_; // in decomposed form!
    Pointer<uint64_t> data_digest_;

    Pointer<uint64_t> data_delta_;

};

struct BatchSelect {
public:

    BatchSelect(const SEALContext::ContextData &context_data, shared_ptr<UniformRandomGenerator> prng) : context_data_(context_data), prng(prng), lhe(context_data, prng), lenc(context_data, prng) {}

    void setup();
    void save_pp(FILE* f) {
        lhe.save_pp(f);
        lenc.save_pp(f);
    }
    void read_pp(FILE* f) {
        lhe.read_pp(f);
        lenc.read_pp(f);
    }

    void enc1(Pointer<uint64_t> &l1);
    void save_st1(FILE* f) {
        lhe.save_st1(f);
    }
    void save_ct1(FILE* f) {
        lhe.save_ct1(f);
        lenc.save_ct1(f);
    }
    void read_st1(FILE* f) {
        lhe.read_st1(f);
    }
    void read_ct1(FILE* f) {
        lhe.read_ct1(f);
        lenc.read_ct1(f);
    }

    void enc2(Pointer<uint64_t> &l2);
    void save_st2(FILE* f) {
        lhe.save_st2(f);
    }
    void read_st2(FILE* f) {
        lhe.read_st2(f);
    }
    void save_ct2(FILE* f) {
        lhe.save_ct2(f);
    }
    void read_ct2(FILE* f) {
        lhe.read_ct2(f);
    }

    void keygen(Pointer<uint64_t> &y);
    void save_sk(FILE* f) {
        lhe.save_sk(f);
    }
    void read_sk(FILE* f) {
        lhe.read_sk(f);
    }

    void dec(Pointer<uint64_t> &y, Pointer<uint64_t> &out);

//private:
    const SEALContext::ContextData &context_data_;
    shared_ptr<UniformRandomGenerator> prng;

    LHE lhe;
    Lenc lenc;

};