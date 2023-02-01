#include <compressed/CBlockTree.h>
#include <pointer_based/BlockTree.h>
#include <unordered_set>
#include <vector>

using Level = BlockTree::Level;

namespace cbt_util {

///
/// @brief Populates the ranks as well as the number of characters inside each block for the first level
/// (beneath the root) of the compressed block tree.
///
/// This corresponds to the fields `CBlockTree::first_level_ranks_` and entry for the first level in
/// `CBlockTree::pop_counts_`.
///
/// @param cbt The compressed block tree whose fields to populate.
/// @param first_level The blocks of the lowest level in the original block tree that has no blocks missing
///
void populate_first_level_ranks(CBlockTree *cbt, const Level &first_level) {
    /// Saves the current cumulative ranks for each character up to the block in the current iteration
    std::unordered_map<int, int> current_prefix_ranks;

    // Iterate through the rank values for the first block to see which characters actually exist in the blocks and
    // create a new empty ranks vector for each of them The new rank vectors are used for the compressed block tree
    for (const auto &[character, _] : first_level[0]->pop_counts_) {
        cbt->first_level_ranks_[character] = new sdsl::int_vector<>(first_level.size());
        cbt->pop_counts_[character].push_back(new sdsl::int_vector<>(first_level.size()));
    }

    // iterate through every block of the lowest level and for each character,
    // calculate the ranks up to this block, as well as the ranks inside this block
    for (int block_index = 0; block_index < first_level.size(); ++block_index) {
        // For every character, save its prefix rank value up to this block into the new map
        // This will take the prefix rank values calculated in the previous iteration,
        // since we don't want to save the rank values *before* this block
        for (const auto &[character, rank] : current_prefix_ranks) {
            (*cbt->first_level_ranks_[character])[block_index] = rank;
        }

        // Lookup the ranks for each character inside the current block
        for (const auto &[character, _] : first_level[block_index]->pop_counts_) {
            // Set the ranks of this character for the current block
            (*cbt->pop_counts_[character][0])[block_index] = first_level[block_index]->pop_counts_[character];
            // Store the cumulative ranks before this block into the temporary map
            // This will be the prefix ranks value for the *next* block (since we only want to count the ranks *before*)
            // each block
            current_prefix_ranks[character] += first_level[block_index]->pop_counts_[character];
        }
    }

    // Compress the rank int vectors
    for (auto &[_, ranks] : cbt->first_level_ranks_) {
        sdsl::util::bit_compress(*ranks);
    }
    for (auto &[_, ranks] : cbt->pop_counts_) {
        sdsl::util::bit_compress(*(ranks[0]));
    }
}

///
/// @brief Populates the pop count values inside each block for the given level.
///
/// This modifies the pop_counts_ field. The rank information in the given level
/// is pushed on top of pop_counts_[character] for each character.
///
/// @param cbt The compressed block tree whose fields to populate.
/// @param level The blocks from the original block tree of the level that is to be processed.
///
void populate_level_pop_counts(CBlockTree *cbt, const Level &level) {
    // Allocate a new int vector for each character to store its ranks
    // We are pushing the new int vector onto the vector containing all int vectors
    // That means that from now on the last element (.back()) for each character is the one we need to worry about for
    // this level
    for (auto &[character, _] : level[0]->pop_counts_) {
        cbt->pop_counts_[character].push_back(new sdsl::int_vector<>(level.size()));
    }

    // Populate the pop count map with the ranks of the characters inside each block of the current level
    for (int i = 0; i < level.size(); ++i) {
        for (auto &[character, rank] : level[i]->pop_counts_) {
            (*cbt->pop_counts_[character].back())[i] = rank;
        }
    }

    // Compress the rank int vectors
    for (auto &[character, ranks] : cbt->pop_counts_) {
        sdsl::util::bit_compress(*ranks.back());
    }
}

/// For every non-internal block
using Offsets   = sdsl::int_vector<>;
using PopCounts = std::unordered_map<int, sdsl::int_vector<> *>;

///
/// @brief Calculate the offsets and first block pop counts for the non-internal blocks on this level.
///
/// For every non-internal block the offsets are the index in the source text from which the block copies its content (=
/// its source). The indices in the int vector skip internal blocks entirely. See also the offsets_ attribute in
/// CBlockTree.
///
/// For every non-internal block the pop first_block_pop_counts are the number of times a character appears in the
/// source of the block, but only in the part which lies in the first block which the block is copying from.
///
/// @param current_level The level to be processed.
/// @param number_of_leaves The number of leaves on this level. Leaves refer to non-internal blocks (so back blocks are
/// included here).
/// @param is_internal_block A bit vector saving for each block whether it is internal or not.
///
std::pair<Offsets *, PopCounts> offsets_and_first_block_pop_counts(const Level      &current_level,
                                                                   const int         number_of_leaves,
                                                                   sdsl::bit_vector &is_internal_block) {
    // For each non-internal of this level, save the total offset (in characters) from the start of the string, from
    // which the content of this block is copied
    auto *offsets = new Offsets(number_of_leaves);
    // For each character, how often does this character appear in the part of the source that lies in first_block
    // See the docs for the respective field in `Block`.
    // The int vectors in here store information for each block of this level that is not an internal block i.e. only
    // for leaf blocks and back blocks
    PopCounts first_block_pop_counts;

    // For every character that exists create a new int vector to store the pop counts
    for (const auto [character, _] : current_level[0]->pop_counts_) {
        first_block_pop_counts[character] = new sdsl::int_vector<>(number_of_leaves);
    }

    // The index of the current block in this level, not counting internal nodes
    int       j            = 0;
    const int block_length = current_level.front()->length();
    for (int i = 0; i < current_level.size(); ++i) {
        // If the current block is not internal (i.e. is a leaf block or a back block)
        if (!is_internal_block[i]) {
            // This is now the j-th non-internal block on this level
            for (const auto [character, first_block_pop_count] : current_level[i]->first_block_pop_counts_) {
                (*first_block_pop_counts[character])[j] = first_block_pop_count;
            }
            // Calculate the offset from the start of the string from which this block copies.
            (*offsets)[j] = current_level[i]->first_block_->level_index_ * block_length + current_level[i]->offset_;
            j++;
        }
    }

    sdsl::util::bit_compress(*offsets);
    return {offsets, first_block_pop_counts};
}

///
/// @brief Fill a level of the the block tree with data used for fast substring queries.
/// This saves the first and last few characters in each block to facilitate substring queries.
///
/// @param cbt The compressed block tree that's under construction
/// @param level The level to process
///
void fill_fast_substring_data(CBlockTree *cbt, const Level &level) {
    const size_t level_size           = level.size();
    const size_t level_block_size     = level[0]->length();
    const size_t saved_data_per_block = std::min(level_block_size, 2 * cbt->prefix_suffix_size_);
    // Create a new int vector to hold this level's symbols
    cbt->prefix_suffix_symbols_.push_back(new sdsl::int_vector<>(level_size * saved_data_per_block, 0));
    auto &symbols = *cbt->prefix_suffix_symbols_.back();
    for (int block_index = 0; block_index < level_size; block_index++) {
        Block           *block  = level[block_index];
        std::string_view prefix = block->prefix_;
        std::string_view suffix = block->suffix_;

        // In this case the prefix and suffix are equal
        if (level_block_size <= cbt->prefix_suffix_size_) {
            for (size_t i = 0; i < saved_data_per_block; ++i) {
                // Get the char from the prefix. If a char does not exist there, fill up the rest with garbage
                // This can only happen for the last block of a level so the space consumed by this is insignificant
                const auto prefix_char =
                    i < prefix.size() && block->start_index_ + i < cbt->input_size_ ? prefix[i] : prefix[0];
                symbols[block_index * saved_data_per_block + i] = cbt->mapping_.at(prefix_char);
            }
            continue;
        }

        // In this case, prefix and suffix overlap
        if (level_block_size < 2 * cbt->prefix_suffix_size_) {
            const size_t overlap = prefix.size() + suffix.size() - saved_data_per_block;
            for (size_t i = 0; i < saved_data_per_block; ++i) {
                // Read characters from the prefix until there are none left
                if (i < prefix.size()) {
                    const char prefix_char =
                        i < prefix.size() && block->start_index_ + i < cbt->input_size_ ? prefix[i] : prefix[0];
                    symbols[block_index * saved_data_per_block + i] = cbt->mapping_.at(prefix_char);
                } else {
                    // From that point on, start reading from the suffix, skipping the characters that overlap
                    const size_t index_in_suffix = i - prefix.size() + overlap;
                    const char   suffix_char =
                        index_in_suffix < suffix.size() && block->start_index_ + i < cbt->input_size_
                              ? suffix[index_in_suffix]
                              : prefix[0];
                    symbols[block_index * saved_data_per_block + i] = cbt->mapping_.at(suffix_char);
                }
            }
            continue;
        }

        for (size_t i = 0; i < saved_data_per_block >> 1; ++i) {
            symbols[block_index * saved_data_per_block + i] = cbt->mapping_[prefix[i]];
        }

        for (size_t i = 0; i < saved_data_per_block >> 1; ++i) {
            // Get the char from the suffix. If a char does not exist there, fill up the rest with garbage
            // This can only happen for the last block of a level so the space consumed by this is insignificant
            const auto suffix_char                                          = i < suffix.size() ? suffix[i] : prefix[0];
            symbols[block_index * saved_data_per_block + prefix.size() + i] = cbt->mapping_.at(suffix_char);
        }
    }
    sdsl::util::bit_compress(symbols);
}

} // namespace cbt_util

CBlockTree::CBlockTree(BlockTree *bt) :
    arity_(bt->arity_),
    root_arity_(bt->root_arity_),
    prefix_suffix_size_(bt->prefix_suffix_size_),
    input_size_(bt->input_.size()),
    rank_select_support_(bt->rank_select_support_) {
    using namespace cbt_util;

    // Get the lowest level in the block tree that has no "gaps"
    // This will de-facto become the first level below the root
    const auto [lowest_complete_level, _] = bt->get_lowest_complete_level();

    // Populate the rank values in the lowest level of the tree that is still complete (i.e. there are no blocks
    // missing). We will use this as the first level beneath the root of the compressed trie.
    populate_first_level_ranks(this, lowest_complete_level);

    first_level_block_size_ = lowest_complete_level[0]->length();
    number_of_levels_       = 0;

    // Continue further down the tree.
    Level current_level = lowest_complete_level;
    Level next_level    = bt->next_level(lowest_complete_level);

    while (!next_level.empty()) {
        // A bit vector for this level that marks whether the specific block is internal (= 1) or not (= 0).
        auto *is_internal_block = new sdsl::bit_vector(current_level.size(), 0);

        int number_of_leaves = 0;
        // Iterate through this level and update each block's level index (this block's index inside the current level),
        // check whether they are internal blocks and count the number of leaves
        for (int i = 0; i < current_level.size(); ++i) {
            current_level[i]->level_index_ = i;

            (*is_internal_block)[i] = !current_level[i]->is_leaf();
            if (current_level[i]->is_leaf()) {
                ++number_of_leaves;
            }
        }

        // Now we populate the pop counts in this level
        populate_level_pop_counts(this, next_level);
        // ...and the offsets as well as the pop counts in first_block_
        auto [offsets, first_block_pop_counts] =
            offsets_and_first_block_pop_counts(current_level, number_of_leaves, *is_internal_block);
        // Add the offsets for this level to the tree
        offsets_.push_back(offsets);

        // Compress the first block pop counts and add them to the tree for each character
        for (auto &[character, pop_counts] : first_block_pop_counts) {
            sdsl::util::bit_compress(*pop_counts);
            (first_block_pop_counts_[character]).push_back(pop_counts);
        }

        // Add the is_internal bitvectors to the tree along with a rank data structure
        is_internal_.push_back(is_internal_block);
        auto *internal_block_ranks = new sdsl::rank_support_v<1>(is_internal_block);
        is_internal_ranks_.push_back(internal_block_ranks);

        // Go to the next level
        current_level = next_level;
        next_level    = bt->next_level(current_level);
        ++number_of_levels_;
    }

    ++number_of_levels_;

    Level last_level = current_level;

    // Calculate the string that is represented by the concatenation of the last level of the tree (all leaves)
    std::string leaf_string;
    for (Block *b : last_level) {
        leaf_string += b->represented_string();
    }

    // Calculate the alphabet
    std::unordered_set<char> alphabet;
    for (char c : leaf_string) {
        alphabet.insert(c);
    }

    // compact the alphabet and create a mapping
    alphabet_   = new sdsl::int_vector<>(alphabet.size());
    int counter = 0;
    for (char c : alphabet) {
        mapping_[c]             = counter;
        (*alphabet_)[counter++] = c;
    }

    sdsl::util::bit_compress(*alphabet_);

    // make the leaf string an int vector
    leaf_string_ = new sdsl::int_vector<>(leaf_string.size());

    for (int i = 0; i < (*leaf_string_).size(); ++i) {
        (*leaf_string_)[i] = mapping_[leaf_string[i]];
    }

    // Generate data for fast substring queries
    if (prefix_suffix_size_ > 0) {
        Level level = lowest_complete_level;
        while (!level.empty()) {
            fill_fast_substring_data(this, level);
            level = bt->next_level(level);
        }
    }

    sdsl::util::bit_compress(*leaf_string_);
}

CBlockTree::CBlockTree(std::istream &in) {
    in.read((char *) &arity_, sizeof(int));
    in.read((char *) &root_arity_, sizeof(int));
    in.read((char *) &first_level_block_size_, sizeof(int));
    in.read((char *) &number_of_levels_, sizeof(int));
    in.read((char *) &rank_select_support_, sizeof(bool));
    in.read((char *) &input_size_, sizeof(size_t));
    in.read((char *) &prefix_suffix_size_, sizeof(size_t));

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        is_internal_.push_back(new sdsl::bit_vector());
        is_internal_[i]->load(in);
    }

    for (sdsl::bit_vector *bv : is_internal_) {
        is_internal_ranks_.push_back(new sdsl::rank_support_v<1>(bv));
    }

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        offsets_.push_back(new sdsl::int_vector<>());
        offsets_[i]->load(in);
    }

    leaf_string_ = new sdsl::int_vector<>();
    leaf_string_->load(in);

    alphabet_ = new sdsl::int_vector<>();
    alphabet_->load(in);

    int c = 0;
    for (int character : (*alphabet_)) {
        mapping_[character] = c++;
    }

    if (prefix_suffix_size_ > 0) {
        for (int i = 0; i < number_of_levels_ - 1; ++i) {
            prefix_suffix_symbols_.push_back(new sdsl::int_vector<>());
            prefix_suffix_symbols_[i]->load(in);
        }
    }

    if (rank_select_support_) {
        for (int character : *alphabet_) {
            first_level_ranks_[character] = new sdsl::int_vector<>();
            first_level_ranks_[character]->load(in);
        }

        for (int character : *alphabet_) {
            for (int i = 0; i < number_of_levels_; ++i) {
                pop_counts_[character].push_back(new sdsl::int_vector<>());
                pop_counts_[character][i]->load(in);
            }
        }

        for (int character : *alphabet_) {
            for (int i = 0; i < number_of_levels_ - 1; ++i) {
                first_block_pop_counts_[character].push_back(new sdsl::int_vector<>());
                first_block_pop_counts_[character][i]->load(in);
            }
        }
    }
}

CBlockTree::~CBlockTree() {

    for (sdsl::bit_vector *bv : is_internal_) {
        delete bv;
    }

    for (sdsl::rank_support_v<1> *rank : is_internal_ranks_) {
        delete rank;
    }

    for (sdsl::int_vector<> *offsets : offsets_) {
        delete offsets;
    }

    for (sdsl::int_vector<> *symbols : prefix_suffix_symbols_) {
        delete symbols;
    }

    delete leaf_string_;
    delete alphabet_;

    for (auto pair : first_level_ranks_) {
        delete pair.second;
    }

    for (const auto &pair : pop_counts_) {
        for (sdsl::int_vector<> *ranks : pair.second) {
            delete ranks;
        }
    }

    for (const auto &pair : first_block_pop_counts_) {
        for (sdsl::int_vector<> *ranks : pair.second) {
            delete ranks;
        }
    }
}

int CBlockTree::access(int i) const {

    int current_block  = i / first_level_block_size_;
    int current_length = first_level_block_size_;
    i -= current_block * first_level_block_size_;
    int level = 0;
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            current_length /= arity_;
            int child_number = i / current_length;
            i -= child_number * current_length;
            current_block = (*is_internal_ranks_[level])(current_block) *arity_ + child_number;
            ++level;
        } else { // Case BackBlock
            int encoded_offset = (*offsets_[level])[current_block - (*is_internal_ranks_[level])(current_block + 1)];
            current_block      = encoded_offset / current_length;
            i += encoded_offset % current_length;
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            }
        }
    }
    return (*alphabet_)[(*leaf_string_)[i + current_block * current_length]];
}

char *CBlockTree::substr_internal(char *buf, size_t idx, size_t len) const {
    size_t index          = idx;
    size_t current_block  = index / first_level_block_size_;
    size_t current_length = first_level_block_size_;
    // Make the index local to the block we are in
    index -= current_block * first_level_block_size_;
    size_t level                      = 0;
    size_t current_prefix_suffix_size = std::min(current_length, (size_t) prefix_suffix_size_);
    size_t saved_data_per_block =
        current_length <= 2 * prefix_suffix_size_ ? current_length : 2 * current_prefix_suffix_size;

    while (level < number_of_levels_ - 1) {
        // If the substring is part of this block's prefix, we can just read it from there
        if (index + len <= current_prefix_suffix_size) {
            // The start index of the substring inside the saves prefix/suffix vector of this level
            const size_t substr_start_index_ = current_block * saved_data_per_block + index;
            for (size_t i = 0; i < len; ++i) {
                const uint16_t unmapped_symbol = (*prefix_suffix_symbols_[level])[substr_start_index_ + i];
                *buf++                         = (*alphabet_)[unmapped_symbol];
            }
            return buf;
        }

        // If the substring is part of this block's suffix, we can just read it from there
        if (index >= current_length - current_prefix_suffix_size && index + len <= current_length) {
            // Make the index local to the suffix
            index -= current_length - current_prefix_suffix_size;
            // We can assume that the block is actually larger than prefix_suffix_size_ at this point and therefore, the
            // prefix and suffix are distinct and saved separately. Because if the block was smaller than
            // prefix_suffix_size_ the prefix and suffix would be identical and only saved once.
            // In that case, we would have already entered the previous if statement and returned
            const size_t suffix_offset      = saved_data_per_block - current_prefix_suffix_size;
            const size_t substr_start_index = current_block * saved_data_per_block + suffix_offset + index;
            for (size_t i = 0; i < len; ++i) {
                const size_t unmapped_symbol = (*prefix_suffix_symbols_[level])[substr_start_index + i];
                *buf++                       = (*alphabet_)[unmapped_symbol];
            }
            return buf;
        }

        // In this case we are crossing a block boundary
        if (index + len > current_length) {
            // Make the index local to the suffix
            index -= current_length - prefix_suffix_size_;
            // If this block is smaller than the prefix_suffix_size then the suffix (= prefix) has no offset
            // Otherwise it is saved after prefix
            const volatile size_t suffix_offset      = saved_data_per_block - prefix_suffix_size_;
            const volatile size_t substr_start_index = current_block * saved_data_per_block + suffix_offset + index;
            // Read the string from this block's suffix and the next block's prefix
            // Since they are contiguous in memory, read them in one go
            for (int i = 0; i < len; ++i) {
                const uint16_t unmapped_symbol = (*prefix_suffix_symbols_[level])[substr_start_index + i];
                *buf++                         = (*alphabet_)[unmapped_symbol];
            }
            return buf;
        }

        // If none of the previous cases apply, the substring is entirely contained in a child block
        // If we are in an internal block, we just go to that child
        // If we are in a back block, we first need to go to its source
        if ((*is_internal_[level])[current_block]) { // Case Internal Block
            // Since each block is split into `arity_` children, this will be the childrens' length
            current_length /= arity_;
            // The number of children we skip by going down in the tree
            int child_number = index / current_length;
            // Skip the children
            index -= child_number * current_length;
            // Since when going down to the children, the back blocks of the current leven are not represented beyond
            // here With the rank query we can ensure that current_block skips the (nonexistent) children of the back
            // blocks
            current_block = (*is_internal_ranks_[level])(current_block) *arity_ + child_number;

            current_prefix_suffix_size = std::min(current_length, current_prefix_suffix_size);
            saved_data_per_block =
                current_length <= 2 * prefix_suffix_size_ ? current_length : 2 * current_prefix_suffix_size;
            level++;
            continue;
        } else { // Case Back Block
            // Find the index of the current block (skipping back blocks)
            auto rank = (*is_internal_ranks_[level])(current_block + 1);
            // Find the position of this block's source in the input text
            int encoded_offset = (*offsets_[level])[current_block - rank];
            current_block      = encoded_offset / current_length;
            index += encoded_offset % current_length;
            // If we happen to move into the next block, we adjust index and current_block
            if (index >= current_length) {
                index -= current_length;
                current_block++;
            }
        }
    }
    // If we made it here without returning, we must be in a leaf
    // So just read from the leaf string
    const size_t start_in_leaf = index + current_block * current_length;
    for (int i = 0; i < len; ++i) {
        *buf++ = (*alphabet_)[(*leaf_string_)[start_in_leaf + i]];
    }
    return buf;
}

char *CBlockTree::substr(char *buf, const size_t index, const size_t len) const {
    // This is the length of the prefix/suffix saved for each block
    const size_t current_prefix_suffix_size =
        first_level_block_size_ < prefix_suffix_size_ ? first_level_block_size_ : prefix_suffix_size_;
    // The number of complete chunks that fit in len
    const size_t complete_chunks = len / current_prefix_suffix_size;
    for (int i = 0; i < complete_chunks; ++i) {
        buf = this->substr_internal(buf, index + i * current_prefix_suffix_size, current_prefix_suffix_size);
    }

    const size_t characters_extracted = complete_chunks * current_prefix_suffix_size;
    // In this case there is a rest of the string we have not extracted yet
    if (characters_extracted < len) {
        buf = substr_internal(buf, index + characters_extracted, len - characters_extracted);
    }
    return buf;
}

int CBlockTree::rank(int c, int i) const {

    auto &ranks        = pop_counts_.at(c);
    auto &second_ranks = first_block_pop_counts_.at(c);

    int current_block  = i / first_level_block_size_;
    int current_length = first_level_block_size_;
    i                  = i - current_block * current_length;
    int level          = 0;

    int r = (*first_level_ranks_.at(c))[current_block];
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            current_length /= arity_;
            int child_number = i / current_length;
            i -= child_number * current_length;

            int firstChild = (*is_internal_ranks_.at(level))(current_block) *arity_;
            for (int child = firstChild; child < firstChild + child_number; ++child) r += (*ranks[level + 1])[child];
            current_block = firstChild + child_number;
            ++level;
        } else { // Case BackBlock
            int index          = current_block - (*is_internal_ranks_[level])(current_block + 1);
            int encoded_offset = (*offsets_.at(level))[index];
            current_block      = encoded_offset / current_length;
            i += encoded_offset % current_length;
            r += (*second_ranks.at(level))[index];
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            } else {
                r -= (*ranks.at(level))[current_block];
            }
        }
    }

    i += current_block * current_length;
    int d = mapping_.at(c);
    for (int j = current_block * current_length; j <= i; ++j) {
        if ((*leaf_string_)[j] == d)
            ++r;
    }

    return r;
}

int CBlockTree::select(int c, int k) const {

    auto &ranks                    = pop_counts_.at(c);
    auto &second_ranks             = first_block_pop_counts_.at(c);
    auto &first_level_prefix_ranks = first_level_ranks_.at(c);

    int current_block = (k - 1) / first_level_block_size_;

    int end_block = first_level_prefix_ranks->size() - 1;
    while (current_block != end_block) {
        int m = current_block + (end_block - current_block) / 2;
        int f = (*first_level_prefix_ranks)[m];
        if (f < k) {
            if (end_block - current_block == 1) {
                if ((*first_level_prefix_ranks)[m + 1] < k) {
                    current_block = m + 1;
                }
                break;
            }
            current_block = m;
        } else {
            end_block = m - 1;
        }
    }

    int current_length = first_level_block_size_;
    int s              = current_block * current_length;
    k -= (*first_level_prefix_ranks)[current_block];
    int level = 0;
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            int firstChild          = (*is_internal_ranks_[level])(current_block) *arity_;
            int child               = firstChild;
            int r                   = (*ranks[level + 1])[child];
            int last_possible_child = (firstChild + arity_ - 1 > (*ranks[level + 1]).size() - 1)
                                          ? (*ranks[level + 1]).size() - 1
                                          : firstChild + arity_ - 1;
            while (child < last_possible_child && k > r) { // Border conditions?
                ++child;
                r += (*ranks[level + 1])[child];
            }
            k -= r - (*ranks[level + 1])[child];
            current_length /= arity_;
            s += (child - firstChild) * current_length;
            current_block = child;
            ++level;
        } else { // Case BackBlock
            int index          = current_block - (*is_internal_ranks_[level])(current_block + 1);
            int encoded_offset = (*offsets_[level])[index];
            current_block      = encoded_offset / current_length;

            k -= (*second_ranks[level])[index];
            s -= encoded_offset % current_length;
            if (k > 0) {
                s += current_length;
                ++current_block;
            } else {
                k += (*ranks[level])[current_block];
            }
        }
    }

    int d = mapping_.at(c);
    for (int j = current_block * current_length;; ++j) {
        if ((*leaf_string_)[j] == d)
            --k;
        if (!k)
            return s + j - current_block * current_length;
    }
}

int CBlockTree::get_partial_size() const {
    int fields = sizeof(int) * 3;

    int leaf_string_size = sdsl::size_in_bytes(*leaf_string_);

    int alphabet_size = sdsl::size_in_bytes(*alphabet_);
    int mapping_size  = sizeof(int) * 256;

    int bt_bv_size = sizeof(void *);
    for (sdsl::bit_vector *bv : is_internal_) {
        bt_bv_size += sdsl::size_in_bytes(*bv);
    }

    int bt_bv_rank_size = sizeof(void *);
    for (sdsl::rank_support_v<1> *bvr : is_internal_ranks_) {
        bt_bv_rank_size += sdsl::size_in_bytes(*bvr);
    }

    int bt_offsets_size = sizeof(void *);
    for (sdsl::int_vector<> *offsets : offsets_) {
        bt_offsets_size += sdsl::size_in_bytes(*offsets);
    }

    int prefix_suffix_size = sizeof(void *);
    for (sdsl::int_vector<> *symbols : prefix_suffix_symbols_) {
        prefix_suffix_size += sdsl::size_in_bytes(*symbols);
    }

    return fields + mapping_size + alphabet_size + bt_bv_size + bt_bv_rank_size + bt_offsets_size + leaf_string_size +
           prefix_suffix_size;
}

int CBlockTree::size() const {

    int bt_ranks_total_size = (pop_counts_.size() + 1) * sizeof(void *);
    for (const auto &pair : pop_counts_) {
        int size = 0;
        for (sdsl::int_vector<> *ranks : pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_ranks_total_size += size;
    }

    int bt_prefix_ranks_first_level_size = 0;
    for (const auto &pair : first_level_ranks_) {
        bt_prefix_ranks_first_level_size += sdsl::size_in_bytes(*(pair.second));
    }

    int bt_second_ranks_total_size = (first_block_pop_counts_.size() + 1) * sizeof(void *);
    for (const auto &pair : first_block_pop_counts_) {
        int size = 0;
        for (sdsl::int_vector<> *ranks : pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_second_ranks_total_size += size;
    }

    int partial_total_size = get_partial_size();
    int rank_size          = bt_second_ranks_total_size + bt_ranks_total_size + bt_prefix_ranks_first_level_size;
    return rank_size + partial_total_size;
}

void CBlockTree::serialize(std::ostream &out) const {

    out.write((char *) &arity_, sizeof(int));
    out.write((char *) &root_arity_, sizeof(int));
    out.write((char *) &first_level_block_size_, sizeof(int));
    out.write((char *) &number_of_levels_, sizeof(int));
    out.write((char *) &rank_select_support_, sizeof(bool));
    out.write((char *) &input_size_, sizeof(size_t));
    out.write((char *) &prefix_suffix_size_, sizeof(size_t));

    for (sdsl::bit_vector *bv : is_internal_) {
        bv->serialize(out);
    }

    for (sdsl::int_vector<> *offsets : offsets_) {
        offsets->serialize(out);
    }

    leaf_string_->serialize(out);

    alphabet_->serialize(out);

    if (prefix_suffix_size_ > 0) {
        for (sdsl::int_vector<> *symbols : prefix_suffix_symbols_) {
            symbols->serialize(out);
        }
    }

    if (rank_select_support_) {
        for (int character : *alphabet_) {
            first_level_ranks_.at(character)->serialize(out);
        }

        for (int character : *alphabet_) {
            for (sdsl::int_vector<> *pop_counts : pop_counts_.at(character)) {
                pop_counts->serialize(out);
            }
        }

        for (int character : *alphabet_) {
            for (sdsl::int_vector<> *pop_counts : first_block_pop_counts_.at(character)) {
                pop_counts->serialize(out);
            }
        }
    }
}
