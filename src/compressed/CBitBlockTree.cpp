#include <algorithm>
#include <compressed/CBitBlockTree.h>
#include <ranges>
#include <unordered_set>
#include <vector>

using Level = BlockTree::Level;

namespace cbbt_util {

///
/// @brief Populates the ranks as well as the number of ones inside each block for the first level
/// (beneath the root) of the compressed block tree.
///
/// This corresponds to the fields `CBitBlockTree::first_level_ranks_` and entry for the first level in
/// `CBitBlockTree::pop_counts_`.
///
/// @param cbbt The compressed block tree whose fields to populate.
/// @param first_level The blocks of the lowest level in the original block tree that has no blocks missing
///
void populate_first_level_ranks(CBitBlockTree *cbbt, const Level &first_level, const int one_symbol) {
    int   current_first_level_ranks = 0;
    auto *first_level_ranks         = new sdsl::int_vector<>(first_level.size());
    auto *first_level_pop_counts    = new sdsl::int_vector<>(first_level.size());

    for (int i = 0; i < first_level.size(); ++i) {
        (*first_level_ranks)[i] = current_first_level_ranks;

        for (auto pair : first_level[i]->pop_counts_) {
            if (pair.first == one_symbol) {
                (*first_level_pop_counts)[i] = first_level[i]->pop_counts_[pair.first];
                current_first_level_ranks    = current_first_level_ranks + first_level[i]->pop_counts_[pair.first];
            }
        }
    }

    sdsl::util::bit_compress(*(first_level_ranks));
    cbbt->first_level_ranks_ = first_level_ranks;

    sdsl::util::bit_compress(*(first_level_pop_counts));
    cbbt->pop_counts_.push_back(first_level_pop_counts);
}

} // namespace cbbt_util

CBitBlockTree::CBitBlockTree(BlockTree *bt, int one_symbol) : arity_(bt->arity_), input_size_(bt->input_.size()) {
    using namespace cbbt_util;
    auto [first_level, _] = bt->get_lowest_complete_level();

    populate_first_level_ranks(this, first_level, one_symbol);

    first_level_block_size_ = first_level[0]->length();
    number_of_levels_       = 0;

    Level current_level = first_level;
    Level next_level    = bt->next_level(first_level);

    while (next_level.size() != 0) {
        // Bit vector for current level that is 1, if the corresponding block is internal
        sdsl::bit_vector *is_internal_block = new sdsl::bit_vector(current_level.size(), 0);

        sdsl::int_vector<> *next_level_pop_counts = new sdsl::int_vector<>(next_level.size());

        int number_of_leaves = 0;
        int current_length   = current_level.front()->length();
        for (int i = 0; i < current_level.size(); ++i) {
            current_level[i]->level_index_ = i;

            (*is_internal_block)[i] = !current_level[i]->is_leaf();
            if (current_level[i]->is_leaf()) {
                ++number_of_leaves;
            }
        }

        // Populate this level's pop counts
        for (int i = 0; i < next_level.size(); ++i) {
            if (next_level[i]->pop_counts_.contains(one_symbol)) {
                (*next_level_pop_counts)[i] = next_level[i]->pop_counts_.at(one_symbol);
            }
        }

        // Populate offsets and pop counts in first_blocks
        sdsl::int_vector<> *offsets                = new sdsl::int_vector<>(number_of_leaves);
        sdsl::int_vector<> *first_block_pop_counts = new sdsl::int_vector<>(number_of_leaves);

        int j = 0;
        for (int i = 0; i < current_level.size(); ++i) {
            if (!(*is_internal_block)[i]) {
                if (current_level[i]->first_block_pop_counts_.contains(one_symbol)) {
                    (*first_block_pop_counts)[j] = current_level[i]->first_block_pop_counts_.at(one_symbol);
                }

                (*offsets)[j++] =
                    current_level[i]->first_block_->level_index_ * current_length + current_level[i]->offset_;
            }
        }

        sdsl::util::bit_compress(*offsets);
        offsets_.push_back(offsets);

        sdsl::util::bit_compress(*(next_level_pop_counts));
        pop_counts_.push_back(next_level_pop_counts);

        sdsl::util::bit_compress(*(first_block_pop_counts));
        first_block_pop_counts_.push_back(first_block_pop_counts);

        is_internal_.push_back(is_internal_block);
        sdsl::rank_support_v<1> *current_level_bv_rank = new sdsl::rank_support_v<1>(is_internal_block);
        is_internal_ranks_.push_back(current_level_bv_rank);

        current_level = next_level;
        next_level    = bt->next_level(current_level);
        ++number_of_levels_;
    }

    ++number_of_levels_;

    // This level consists only of leaves
    Level last_level = current_level;

    std::string leaf_string = "";
    for (Block *b : last_level) {
        leaf_string += b->represented_string();
    }

    // Create the leaf bit vector from the concatenated string of all leaves
    leaf_bv_ = new sdsl::bit_vector(leaf_string.size());

    for (int i = 0; i < (*leaf_bv_).size(); ++i) {
        (*leaf_bv_)[i] = leaf_string[i] == one_symbol;
    }
}

CBitBlockTree::CBitBlockTree(std::istream &in) {
    in.read((char *) &arity_, sizeof(int));
    in.read((char *) &first_level_block_size_, sizeof(int));
    in.read((char *) &number_of_levels_, sizeof(int));
    in.read((char *) &input_size_, sizeof(size_t));

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

    leaf_bv_ = new sdsl::bit_vector();
    (*leaf_bv_).load(in);

    first_level_ranks_ = new sdsl::int_vector<>();
    (*first_level_ranks_).load(in);

    for (int i = 0; i < number_of_levels_; ++i) {
        pop_counts_.push_back(new sdsl::int_vector<>());
        (*pop_counts_[i]).load(in);
    }

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        first_block_pop_counts_.push_back(new sdsl::int_vector<>());
        (*first_block_pop_counts_[i]).load(in);
    }
}

CBitBlockTree::~CBitBlockTree() {

    for (sdsl::bit_vector *bv : is_internal_) {
        delete bv;
    }

    for (sdsl::rank_support_v<1> *rank : is_internal_ranks_) {
        delete rank;
    }

    for (sdsl::int_vector<> *offsets : offsets_) {
        delete offsets;
    }

    delete leaf_bv_;

    if (first_level_ranks_)
        delete first_level_ranks_;

    for (auto ranks : pop_counts_) {
        delete ranks;
    }

    for (auto ranks : first_block_pop_counts_) {
        delete ranks;
    }
}

int CBitBlockTree::access(int i) const {

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
    return (*leaf_bv_)[i + current_block * current_length];
}

int CBitBlockTree::rank_1(int i) const {

    auto &ranks        = pop_counts_;
    auto &second_ranks = first_block_pop_counts_;

    int current_block  = i / first_level_block_size_;
    int current_length = first_level_block_size_;
    i                  = i - current_block * current_length;
    int level          = 0;

    int r = (*first_level_ranks_)[current_block];
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

    auto    &leaf_bv    = *leaf_bv_;
    int      chunk      = (current_block * current_length) / 64;
    uint64_t chunk_info = *(leaf_bv.m_data + chunk);

    i += current_block * current_length;
    for (int j = current_block * current_length; j <= i; ++j) {
        int value = (chunk_info >> (j % 64)) % 2;
        if (value == 1)
            ++r;
        if ((j + 1) % 64 == 0) {
            ++chunk;
            chunk_info = *(leaf_bv.m_data + chunk);
        }
    }

    return r;
}

int CBitBlockTree::rank_0(int i) const { return i - rank_1(i) + 1; }

int CBitBlockTree::select_1(int k) const {

    auto &ranks                    = pop_counts_;
    auto &second_ranks             = first_block_pop_counts_;
    auto &first_level_prefix_ranks = first_level_ranks_;

    int current_block = (k - 1) / first_level_block_size_;

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
            while (child < last_possible_child && k > r) {
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

    auto    &leaf_bv    = *leaf_bv_;
    int      chunk      = (current_block * current_length) / 64;
    uint64_t chunk_info = *(leaf_bv.m_data + chunk);

    for (int j = current_block * current_length;; ++j) {
        int value = (chunk_info >> (j % 64)) % 2;
        if (value == 1)
            --k;
        if (!k)
            return s + j - current_block * current_length;
        if ((j + 1) % 64 == 0) {
            ++chunk;
            chunk_info = *(leaf_bv.m_data + chunk);
        }
    }

    return -1;
}

int CBitBlockTree::select_0(int k) const {

    auto &ranks                    = pop_counts_;
    auto &second_ranks             = first_block_pop_counts_;
    auto &first_level_prefix_ranks = first_level_ranks_;

    int current_block = (k - 1) / first_level_block_size_;

    int end_block = (*first_level_prefix_ranks).size() - 1;
    while (current_block != end_block) {
        int m = current_block + (end_block - current_block) / 2;
        int f = first_level_block_size_ * m - (*first_level_prefix_ranks)[m];
        if (f < k) {
            if (end_block - current_block == 1) {
                if ((first_level_block_size_ * (m + 1) - (*first_level_prefix_ranks)[m + 1]) < k) {
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
    k -= s - (*first_level_prefix_ranks)[current_block];
    int level = 0;
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            int firstChild          = (*is_internal_ranks_[level])(current_block) *arity_;
            int child               = firstChild;
            int child_length        = current_length / arity_;
            int r                   = child_length - (*ranks[level + 1])[child];
            int last_possible_child = (firstChild + arity_ - 1 > (*ranks[level + 1]).size() - 1)
                                          ? (*ranks[level + 1]).size() - 1
                                          : firstChild + arity_ - 1;
            while (child < last_possible_child && k > r) {
                ++child;
                r += child_length - (*ranks[level + 1])[child];
            }
            k -= r - (child_length - (*ranks[level + 1])[child]);
            current_length = child_length;
            s += (child - firstChild) * current_length;
            current_block = child;
            ++level;
        } else { // Case BackBlock
            int index          = current_block - (*is_internal_ranks_[level])(current_block + 1);
            int encoded_offset = (*offsets_[level])[index];
            current_block      = encoded_offset / current_length;

            k -= (current_length - encoded_offset % current_length) - (*second_ranks[level])[index];
            s -= encoded_offset % current_length;
            if (k > 0) {
                s += current_length;
                ++current_block;
            } else {
                k += current_length - (*ranks[level])[current_block];
            }
        }
    }

    auto    &leaf_bv    = *leaf_bv_;
    int      chunk      = (current_block * current_length) / 64;
    uint64_t chunk_info = *(leaf_bv.m_data + chunk);

    for (int j = current_block * current_length;; ++j) {
        int value = (chunk_info >> (j % 64)) % 2;
        if (value == 0)
            --k;
        if (!k)
            return s + j - current_block * current_length;
        if ((j + 1) % 64 == 0) {
            ++chunk;
            chunk_info = *(leaf_bv.m_data + chunk);
        }
    }

    return -1;
}

int CBitBlockTree::get_partial_size() const {

    int fields = sizeof(int) * 3;

    int leaf_bv_size = sdsl::size_in_bytes(*leaf_bv_);

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

    return fields + bt_bv_size + bt_bv_rank_size + bt_offsets_size + leaf_bv_size;
}

int CBitBlockTree::size() const {

    int bt_ranks_total_size = sizeof(void *);
    for (auto ranks : pop_counts_) {
        bt_ranks_total_size += sdsl::size_in_bytes(*ranks);
    }

    int bt_prefix_ranks_first_level_size = sdsl::size_in_bytes(*(first_level_ranks_));

    int bt_second_ranks_total_size = sizeof(void *);
    for (auto ranks : first_block_pop_counts_) {
        bt_second_ranks_total_size += sdsl::size_in_bytes(*ranks);
    }

    int partial_total_size = get_partial_size();
    int rank_size          = bt_second_ranks_total_size + bt_ranks_total_size + bt_prefix_ranks_first_level_size;
    return rank_size + partial_total_size;
}

void CBitBlockTree::serialize(std::ostream &out) const {

    out.write((char *) &arity_, sizeof(int));
    out.write((char *) &first_level_block_size_, sizeof(int));
    out.write((char *) &number_of_levels_, sizeof(int));
    out.write((char *) &input_size_, sizeof(size_t));

    for (sdsl::bit_vector *bv : is_internal_) {
        (*bv).serialize(out);
    }

    for (sdsl::int_vector<> *offsets : offsets_) {
        (*offsets).serialize(out);
    }

    (*leaf_bv_).serialize(out);

    (*first_level_ranks_).serialize(out);

    for (sdsl::int_vector<> *ranks : pop_counts_) {
        (*ranks).serialize(out);
    }

    for (sdsl::int_vector<> *second_ranks : first_block_pop_counts_) {
        (*second_ranks).serialize(out);
    }
}
