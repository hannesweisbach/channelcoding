channelcoding =============

My stuff for channel coding lecture

I originally started coding this because I kept making arithmetic errors when
doing the exercises for my channel coding lecture.  I then re-used the code for
simulating iterative soft-decoding algorithms on BCH parity check matrices.

In the gf/ subfolder, Galois field arithmetic is implemented.

The 'gf' class implements Galois field arithmetic for GF(2^q). The class is
templated for the power q and the modular polynomial used to construct the
field.  The modular polynomial has defaults for a q of up to 8.

The Galois field element is implemented as inner class of the gf class.  In
conjuction with the template parameters for the gf class, this ensures
type-safety. Binary arithmetic operations can only be performed on elements of
the same field (unless explicit casting, i.e. static_cast<>(), is used).

Alternatively, the field can be modelled as inner class of the field element.
The only reason I need the field explicitly, is to iterate over all elements of
the field. Therefore, the gf class provides [c|r]begin()/end() member
functions. I decided against overloading pre- and/or post-increment
operator++/(int) for the purpose of iterating over the field, because there
would be a semantic disparity. On the one hand, one expects operator++ to be
semantically equivalent to adding 1:

using Element = gf::gf<5>::element_t;

Element a(0); a++;

assert(a == Element(0) + Element(1));

However, in Galois field arithmetic adding 1 repeatedly, gets you nowhere:

using Element = gf::gf<5>::element_t;

Element a(0); a++; a++; a++; a++;

assert(a == Element(0));

To iterate over all elements of a field, one would have to increment the power
of a primitive (generating) element.  Doing this has two drawbacks: 1) It does
not generate the element 0 straight away.  2) The semantic of operator++ has
changed from 'adding 1' to 'multiply with the primitive element of the field'.

Therefore, I refrained from overloading operator++/(int) for this purpose.

The backing storage of an element of a field is chosen depending on the size of
the field.  The fields from GF(2) to GF(2^8) use uint8_t as backing storage,
fields up to GF(2^16) use uint16_t, and so forth.

A Galoid Field element overloads the usual arithmetic operators +, -, *, /, as
well as += and *=.  Also overloads for the relational operators <, ==, and !=
are provided. Providing overloads for the missing operators is left as TODO.
Additional the output stream operator as well as a number of conversion
operators are provided.

Also, it is desireable to be able to create the values 0 and 1 from integer
literals such that generic algorithms can use a template type parameter T to
create the like so: T(0) or T(1).

Polynomials

The polynomial class is templated over the coefficient type. It is intended to
be as general as possible. For this end, algorithms like finding roots are
implemented as free standing functions.

The usual arithmetic operators for polynomial arithmetic are provided as well
as the relational operators == and !=. Operator() is overloaded for evaluating
the polynomial at a given point.

The gcd and lcm are implemented as member functions but should be implemented
as algorithms, with the same reasoning as for finding zeroes.

The member function n(size_t k) is provided to construct the polynomial x^k. A
more intuitive way to construct 'x^k + e'-literals is desirable.

Polynomials derive from std::vector and thus are dynamically sized to be more
generic. For the sole purpose of providing polynomial arithmetic for cyclic
codes polynomials with fixed maximum degree can be used. This would make it
possible to construct polynomials at compile-time and possibly some performance
could be gained by avoiding dynamic memory allocation at the cost of
overprovisioned memory.

Both, Galois field arithmetic as well as polynomial arithmetic could use some
test suites. A specialization for polynomials over GF(2) would also be nice.

Codes

The implementation is focused on cyclic codes; primarily BCH and RS codes.
Encoding and decoding is implemented as generic algorithm for the
multiplication and division methods. The requirement for the Polynomial
template parameter is having a member degree() for determining the degree of
the polynomial, as well as a static member function n(size_t) to construct the
polynomial x^k.

The cyclic class (which should be named bch_base) is templated for the power q
of the extension field, the error correction capability, the correction
algorithm, the length of the code, the encoding method, and the algorithm to
determine the error values.

The correction capability can be specified in terms of correctable errors or
minimum distance. Specifying the length of the code allows for code shortening.

Primitive_bch and rs derive from the cyclic class. They provide code-specific
ways of determining the error values. For binary BCH codes the error values are
always 1, for RS codes the naive algorithm consisting of solving an equation
system of error positions and syndromes is provided. Implementing the more
efficient Forney algorithm is left as a TODO.

The encode() member function takes an InputSequence as parameter. The
requirements on InputSequence are a size() member function, overloads for
std::cbegin()/std::cend() and a value_type type declaration.  Alternatively, an
input iterator pair for begin and end could be used. This is not implemented.
The channel code word is return by means of an output iterator. Objects of type
InputSequence::value_type are assigned to the output iterator.
The reason for using an output iterator is to facilitate zero-copy encoding,
when the channel code word has to be written in a special buffer for
transmission.

The decode() member function takes an InputSequence and an optional
std::vector<unsigned> of erased positions. The decoded source code word is
return via a std::vector<T>, where T is template parameter which can be
specified by the caller but is defaulted to InputSequence::value_type.

To make simulations less computationally expensive by avoiding
encoding/decoding a member function correct() is exposed, only providing
correction but no decoding. Similar to decode(), the parameters for correct()
are an InputSequence and a std::vector of erasures. The corrected sequence is
returned as std::vector<T>, but here T is defaulted to the underlying storage
type of the Galois field elements (which is usually uint8_t).

The reason for these diverging and incnosistent interfaces are a result of me
experimenting with different possiblities of these interfaces and finding an
acceptable trade-off between generality of the algorithm and performance of the
implementation.  Although this code is intended for demonstration purposes,
simulation, and benchmarking the performance of the implementation is of
interest, since a high number of simulations with long codes can amount to
quite some computing time.
This implementation is not production quality code, and more importantly, a
specific use case requires only a specific code, which can be implemented more
efficiently than this general approach.

Additionally, the primitive_bch class provides an overload for the correct()
member function to implement erasure decoding with the
Peterson-Gorenstein-Zierler algorithm.

For finding error values the Berlekamp-Massey algorithm and the Euklid
algorithm are also implemented. By constructing a parity check matrix H for the
BCH codes, iterative (soft-decision) Min-Sum decoding can also be used.

Different modification for Min-Sum decoding are implemted, namely Min-Sum,
offset Min-Sum, normalized Min-Sum, two versions of self-correcting Min-Sum and
2D-normalized Min-Sum.
The Min-Sum algorithm can be modified by template parameters.

For implementing the Min-Sum algorithm a matrix class is used. It is only used
as storage and the only requirement are the member functions at(size_t) and
operator[] accordingly.

Contrary to LDPC codes, the matrix class is not implemented as sparse matrix,
since parity check matrices of BCH codes are not sparse.

The uncoded class implements a no-op channel code to provide a reference.

Simulation

To evaluate code performance channel simulations are done to measure the word
error rate (WER) at different Eb/N0 conditions.  The awgn_simulation class
provides simulation of an AWGN channel, wheras the bitflip_simulation class
provides exhaustive test from 0 up to specified number of bitflips. For each
number of bitflips every possible received word is generated and the correction
algorithm is run.
The implementation of simulation environments of a BSC and a BEC is left as a
TODO.

To avoid the necessity for all code classes to derive from a common base to be
used polynmorphically, the decoder class is used to hide virtual inheritance
between unrelated types. This technique is also used in the Adobe poly<>
library and is called concepts-based polymorphism (Sean Parent).

The uncoded program simulates uncoded, BPSK-modulated communication over an
AWGN channel.

Misc

The exercises program implements solutions for the tasks posed in the lecture.
The tasks are available via
http://www.inf.tu-dresden.de/content/institutes/sya/dud/lectures/2013wintersemester/Kanalkodierung/aufg_samml.pdf

The tasks for admittance to the exam are available via 
http://www.inf.tu-dresden.de/content/institutes/sya/dud/lectures/2014wintersemester/Kanalkodierung/leistungsnachweiseWS1415.pdf

The reuslt is available in the iterative_soft_decoding_of_bch_codes.pdf in this repository.

This implementation requires a C++14 (N3797) conformant compiler and standard
library.  The only conformant library as of this writing is libc++. For OSX
Mountain Lion and later this works out-of-the-box. For Linux libc++ has to be
installed and used.
