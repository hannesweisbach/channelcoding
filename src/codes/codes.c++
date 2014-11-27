#include "codes.h"

decoding_failure::~decoding_failure() = default;
decoding_failure::decoding_failure(const decoding_failure &) = default;
decoding_failure::decoding_failure(decoding_failure &&) = default;
decoding_failure &decoding_failure::operator=(const decoding_failure &) =
    default;
decoding_failure &decoding_failure::operator=(decoding_failure &&) = default;

