#include <compressed/CBlockTree.h>
#include <pointer_based/BlockTree.h>
#include <unordered_set>
#include <vector>

using Level = BlockTree::Level;

///
/// @brief Gets the lowest level of the given block tree that is not fully comprised of internal blocks.
///
/// In other words, this returns the lowest level of the block tree that is still *complete* in the sense that there are
/// no gaps in the level where blocks are missing.
///
/// @param bt The uncompressed pointer-based block tree
/// @return A vector of pointers pointing to the blocks of the lowest complete level of the given block tree in order
/// from left to right.
///
auto get_lowest_complete_level(const BlockTree *bt) -> Level {
    Level current_level = {bt->root_block_};
    while (true) {
        for (Block *b : current_level) {
            if (b->is_leaf()) {
                return current_level;
            }
        }
        current_level = bt->next_level(current_level);
    }
}

///
/// @brief Populates the ranks as well as the number of characters inside each block for the lowest complete level of
/// the compressed block tree.
///
/// This corresponds to the fields `CBlockTree::lowest_complete_level_ranks_` and entry for the lowest complete level in
/// `CBlockTree::pop_counts_`.
///
/// @param cbt The compressed block tree whose fields to populate.
/// @param lowest_complete_level The blocks of the lowest complete level in the original block tree in order from
/// left to right.
///
void populate_lowest_complete_level_ranks(CBlockTree *cbt, const Level &lowest_complete_level) {
    /// Saves the current cumulative ranks for each character up to the block in the current iteration
    std::unordered_map<int, int> current_prefix_ranks;

    // Iterate through the rank values for the first block to see which characters actually exist in the blocks and
    // create a new empty ranks vector for each of them The new rank vectors are used for the compressed block tree
    for (const auto &[character, _] : lowest_complete_level[0]->pop_counts_) {
        cbt->lowest_complete_level_ranks_[character] = new sdsl::int_vector<>(lowest_complete_level.size());
        cbt->pop_counts_[character].push_back(new sdsl::int_vector<>(lowest_complete_level.size()));
    }

    // iterate through every block of the lowest level and for each character,
    // calculate the ranks up to this block, as well as the ranks inside this block
    for (int block_index = 0; block_index < lowest_complete_level.size(); ++block_index) {
        // For every character, save its prefix rank value up to this block into the new map
        // This will take the prefix rank values calculated in the previous iteration,
        // since we don't want to save the rank values *before* this block
        for (const auto &[character, rank] : current_prefix_ranks) {
            (*cbt->lowest_complete_level_ranks_[character])[block_index] = rank;
        }

        // Lookup the ranks for each character inside the current block
        for (const auto &[character, _] : lowest_complete_level[block_index]->pop_counts_) {
            // Set the ranks of this character for the current block
            (*cbt->pop_counts_[character][0])[block_index] = lowest_complete_level[block_index]->pop_counts_[character];
            // Store the cumulative ranks before this block into the temporary map
            // This will be the prefix ranks value for the *next* block (since we only want to cound the ranks *before*)
            // each block
            current_prefix_ranks[character] += lowest_complete_level[block_index]->pop_counts_[character];
        }
    }

    // Compress the rank intvectors
    for (auto &[_, ranks] : cbt->lowest_complete_level_ranks_) {
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
    // We are pushing the new int vector onto the vector containing all intvectors
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
                                                                   sdsl::bit_vector *is_internal_block) {
    // For each non-internal of this level, save the total offset (in characters) from the start of the string, from
    // which the content of this block is copied
    Offsets *offsets = new sdsl::int_vector<>(number_of_leaves);
    // For each character, how often does this character appear in the part of the source that lies in first_block
    // See the docs for the respective field in `Block`.
    // The intvectors in here store information for each block of this level that is not an internal block i.e. only
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
        if (!(*is_internal_block)[i]) {
            // This is now the j-th non-internal block on this level
            for (const auto [character, first_block_pop_count] : current_level[i]->pop_counts_in_first_block_) {
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

CBlockTree::CBlockTree(BlockTree *bt) : arity_(bt->arity_), root_arity_(bt->root_arity_) {
    Level lowest_complete_level = get_lowest_complete_level(bt);

    // Populate the rank values in the lowest level of the tree that is still complete (i.e. there are no blocks
    // missing)
    populate_lowest_complete_level_ranks(this, lowest_complete_level);

    lowest_complete_level_block_size_ = lowest_complete_level[0]->length();
    number_of_levels_                 = 0;

    // Continue further down the tree.
    Level current_level = lowest_complete_level;
    Level next_level    = bt->next_level(lowest_complete_level);

    while (next_level.size() != 0) {
        // A bit vector for this level that marks whether the specific block is internal (= 1) or not (= 0).
        sdsl::bit_vector *is_internal_block = new sdsl::bit_vector(current_level.size(), 0);

        int number_of_leaves = 0;
        // Iterate through this level and update each block's level index (this block's index inside the current level),
        // check whether they are internal blocks and count the number of leaves
        for (int i = 0; i < current_level.size(); ++i) {
            current_level[i]->level_index_ = i;

            if (current_level[i]->is_leaf()) {
                (*is_internal_block)[i] = 0;
                ++number_of_leaves;
            } else {
                (*is_internal_block)[i] = 1;
            }
        }

        // Now we populate the pop counts in this level
        populate_level_pop_counts(this, next_level);
        // ...and the offsets as well as the pop counts in first_block_
        auto [offsets, first_block_pop_counts] =
            offsets_and_first_block_pop_counts(current_level, number_of_leaves, is_internal_block);
        // Add the offsets for this level to the tree
        offsets_.push_back(offsets);

        // Compress the first block pop counts and add them to the tree for each characters
        for (auto &[character, pop_counts] : first_block_pop_counts) {
            sdsl::util::bit_compress(*pop_counts);
            (first_block_pop_counts_[character]).push_back(pop_counts);
        }

        // Add the is_internal bitvectors to the tree along with a rank data structure
        is_internal_.push_back(is_internal_block);
        sdsl::rank_support_v<1> *internal_block_ranks = new sdsl::rank_support_v<1>(is_internal_block);
        is_internal_ranks_.push_back(internal_block_ranks);

        // Go to the next level
        current_level = next_level;
        next_level    = bt->next_level(current_level);
        ++number_of_levels_;
    }

    ++number_of_levels_;

    Level last_level = current_level;

    // Calculate the string that is represented by the concatenation of the last level of the tree (all leaves)
    std::string leaf_string = "";
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

    sdsl::util::bit_compress(*leaf_string_);
}

CBlockTree::CBlockTree(std::istream &in) {
    in.read((char *) &arity_, sizeof(int));
    in.read((char *) &root_arity_, sizeof(int));
    in.read((char *) &lowest_complete_level_block_size_, sizeof(int));
    in.read((char *) &number_of_levels_, sizeof(int));
    in.read((char *) &rank_select_support_, sizeof(bool));

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        is_internal_.push_back(new sdsl::bit_vector());
        (*is_internal_[i]).load(in);
    }

    for (sdsl::bit_vector *bv : is_internal_) {
        is_internal_ranks_.push_back(new sdsl::rank_support_v<1>(bv));
    }

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        offsets_.push_back(new sdsl::int_vector<>());
        (*offsets_[i]).load(in);
    }

    leaf_string_ = new sdsl::int_vector<>();
    (*leaf_string_).load(in);

    alphabet_ = new sdsl::int_vector<>();
    (*alphabet_).load(in);

    int c = 0;
    for (int character : (*alphabet_)) {
        mapping_[character] = c++;
    }

    if (rank_select_support_) {
        for (int character : (*alphabet_)) {
            lowest_complete_level_ranks_[character] = new sdsl::int_vector<>();
            (*lowest_complete_level_ranks_[character]).load(in);
        }

        for (int character : (*alphabet_)) {
            for (int i = 0; i < number_of_levels_; ++i) {
                pop_counts_[character].push_back(new sdsl::int_vector<>());
                (*pop_counts_[character][i]).load(in);
            }
        }

        for (int character : (*alphabet_)) {
            for (int i = 0; i < number_of_levels_ - 1; ++i) {
                first_block_pop_counts_[character].push_back(new sdsl::int_vector<>());
                (*first_block_pop_counts_[character][i]).load(in);
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

    delete leaf_string_;
    delete alphabet_;

    if (rank_select_support_) {
        for (auto pair : lowest_complete_level_ranks_) {
            delete pair.second;
        }

        for (auto pair : pop_counts_) {
            for (sdsl::int_vector<> *ranks : pair.second) {
                delete ranks;
            }
        }

        for (auto pair : first_block_pop_counts_) {
            for (sdsl::int_vector<> *ranks : pair.second) {
                delete ranks;
            }
        }
    }
}

int CBlockTree::access(int i) {

    int current_block  = i / lowest_complete_level_block_size_;
    int current_length = lowest_complete_level_block_size_;
    i -= current_block * lowest_complete_level_block_size_;
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

int CBlockTree::rank(int c, int i) {

    auto &ranks        = pop_counts_[c];
    auto &second_ranks = first_block_pop_counts_[c];

    int current_block  = i / lowest_complete_level_block_size_;
    int current_length = lowest_complete_level_block_size_;
    i                  = i - current_block * current_length;
    int level          = 0;

    int r = (*lowest_complete_level_ranks_[c])[current_block];
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            current_length /= arity_;
            int child_number = i / current_length;
            i -= child_number * current_length;

            int firstChild = (*is_internal_ranks_[level])(current_block) *arity_;
            for (int child = firstChild; child < firstChild + child_number; ++child) r += (*ranks[level + 1])[child];
            current_block = firstChild + child_number;
            ++level;
        } else { // Case BackBlock
            int index          = current_block - (*is_internal_ranks_[level])(current_block + 1);
            int encoded_offset = (*offsets_[level])[index];
            current_block      = encoded_offset / current_length;
            i += encoded_offset % current_length;
            r += (*second_ranks[level])[index];
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            } else {
                r -= (*ranks[level])[current_block];
            }
        }
    }

    i += current_block * current_length;
    int d = mapping_[c];
    for (int j = current_block * current_length; j <= i; ++j) {
        if ((*leaf_string_)[j] == d)
            ++r;
    }

    return r;
}

int CBlockTree::select(int c, int k) {

    auto &ranks                    = pop_counts_[c];
    auto &second_ranks             = first_block_pop_counts_[c];
    auto &first_level_prefix_ranks = lowest_complete_level_ranks_[c];

    int current_block = (k - 1) / lowest_complete_level_block_size_;

    int end_block = (*first_level_prefix_ranks).size() - 1;
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

    int current_length = lowest_complete_level_block_size_;
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

    int d = mapping_[c];
    for (int j = current_block * current_length;; ++j) {
        if ((*leaf_string_)[j] == d)
            --k;
        if (!k)
            return s + j - current_block * current_length;
    }

    return -1;
}

int CBlockTree::get_partial_size() {
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

    return fields + mapping_size + alphabet_size + bt_bv_size + bt_bv_rank_size + bt_offsets_size + leaf_string_size;
}

int CBlockTree::size() {

    int bt_ranks_total_size = (pop_counts_.size() + 1) * sizeof(void *);
    for (auto pair : pop_counts_) {
        int size = 0;
        for (sdsl::int_vector<> *ranks : pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_ranks_total_size += size;
    }

    int bt_prefix_ranks_first_level_size = 0;
    for (auto pair : lowest_complete_level_ranks_) {
        bt_prefix_ranks_first_level_size += sdsl::size_in_bytes(*(pair.second));
    }

    int bt_second_ranks_total_size = (first_block_pop_counts_.size() + 1) * sizeof(void *);
    for (auto pair : first_block_pop_counts_) {
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

void CBlockTree::serialize(std::ostream &out) {

    out.write((char *) &arity_, sizeof(int));
    out.write((char *) &root_arity_, sizeof(int));
    out.write((char *) &lowest_complete_level_block_size_, sizeof(int));
    out.write((char *) &number_of_levels_, sizeof(int));
    out.write((char *) &rank_select_support_, sizeof(bool));

    for (sdsl::bit_vector *bv : is_internal_) {
        (*bv).serialize(out);
    }

    for (sdsl::int_vector<> *offsets : offsets_) {
        (*offsets).serialize(out);
    }

    (*leaf_string_).serialize(out);

    (*alphabet_).serialize(out);

    if (rank_select_support_) {
        for (int character : (*alphabet_)) {
            (*lowest_complete_level_ranks_[character]).serialize(out);
        }

        for (int character : (*alphabet_)) {
            for (sdsl::int_vector<> *ranks : pop_counts_[character]) {
                (*ranks).serialize(out);
            }
        }

        for (int character : (*alphabet_)) {
            for (sdsl::int_vector<> *second_ranks : first_block_pop_counts_[character]) {
                (*second_ranks).serialize(out);
            }
        }
    }
}
