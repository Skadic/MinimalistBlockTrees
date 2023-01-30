#ifndef BLOCKTREE_RABINKARP_H
#define BLOCKTREE_RABINKARP_H

#include <string>

///
/// @brief A structure managing Rabin-Karp fingerprints of a string.
///
/// This provides hash values of a text window that can easily be "shifted along" the text using the `next()` method
/// and retrieved using the `hash` method.
///
class RabinKarp {

    /// The value that the hash values are limited to (exclusive).
    /// In practice, this should be a large prime number.
    const uint64_t range_;

    /// The current index in the text for which the current hash is computed
    uint64_t current_index_;
    /// TODO Find out what this is
    uint64_t rm_;

  public:
    /// The alphabet size
    const uint64_t sigma_;
    /// The current hash value
    uint64_t hash_;
    /// The window size for which the hashes are computed
    const uint64_t     window_size_;
    const std::string &s_;

    ///
    /// @brief Creates a new Rabin Karp fingerprinting data structure.
    ///
    /// @param s The input text.
    /// @param start_index The start index in the text from which the computation should start
    /// @param window_size The window_size for which fingerprints should be computed.
    /// @param range The value that the fingerprints are limited to. This should be a large prime number.
    /// @param sigma The alphabet size
    ///
    RabinKarp(const std::string &s, const int start_index, const int window_size, const int range, const int sigma = 257);

    ///
    /// @brief Get the current hash value.
    ///
    /// @return The hash value of the string at index current_index_ with length window_size_.
    ///
    uint64_t hash() const;

    ///
    /// @brief Shifts the index to the left by one and calculates that string's hash value.
    ///
    /// It can then be retrieved using the `hash()` method.
    ///
    void next();
};

#endif // BLOCKTREE_RABINKARP_H
