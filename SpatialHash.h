/*
MIT License

Copyright (c) 2020 Pawel Böning

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SHASH_H
#define SHASH_H

#include <vector>
#include <cmath>

namespace shash
{
    using HashValue = uint32_t;

    using HashFuncPtr = HashValue (*)(void *buf);
    using ReduceFuncPtr = HashValue (*)(HashValue hash, uint32_t range);

    namespace hashing
    {
        struct Murmur
        {
            static const uint32_t SEED = 15953071;

            inline static uint32_t rotl32(uint32_t x, uint8_t r)
            {
                return (x << r) | (x >> (32 - r));
            }

            inline static uint32_t fmix(uint32_t h)
            {
                h ^= h >> 16;
                h *= 0x85ebca6b;
                h ^= h >> 13;
                h *= 0xc2b2ae35;
                h ^= h >> 16;
                return h;
            }

            inline static HashValue hash(void *buf)
            {
                int length = 4;
                uint32_t *b = (uint32_t *)buf;

                uint32_t h1 = SEED;
                uint32_t k1 = 0;

                for (int i = 0; i < length; i++)
                {
                    k1 = (uint32_t)b[i];

                    k1 *= 0xcc9e2d51;
                    k1 = rotl32(k1, 15);
                    k1 *= 0x1b873593;

                    h1 ^= k1;
                    h1 = rotl32(h1, 13);
                    h1 = h1 * 5 + 0xe6546b64;
                }

                h1 ^= (uint32_t)(length * 4);
                h1 = fmix(h1);

                return h1;
            }
        };

        struct xxHash
        {
            static const uint32_t PRIME32_1 = 2654435761U;
            static const uint32_t PRIME32_2 = 2246822519U;
            static const uint32_t PRIME32_3 = 3266489917U;
            static const uint32_t PRIME32_4 = 668265263U;
            static const uint32_t PRIME32_5 = 374761393U;
            static const uint32_t SEED = 15953071;

            inline static uint32_t rotate_left(uint32_t value, int count)
            {
                return (value << count) | (value >> (32 - count));
            }

            inline static uint32_t sub_hash(uint32_t value, uint32_t read_value)
            {
                value += read_value * PRIME32_2;
                value = rotate_left(value, 13);
                value *= PRIME32_1;
                return value;
            }

            inline static HashValue hash(void *buf)
            {
                uint32_t *b = (uint32_t *)buf;
                int len = 4;

                uint32_t h32;
                int index = 0;

                if (len >= 4)
                {
                    int limit = len - 4;
                    uint32_t v1 = SEED + PRIME32_1 + PRIME32_2;
                    uint32_t v2 = SEED + PRIME32_2;
                    uint32_t v3 = SEED + 0;
                    uint32_t v4 = SEED - PRIME32_1;

                    do
                    {
                        v1 = sub_hash(v1, b[index]);
                        index++;
                        v2 = sub_hash(v2, b[index]);
                        index++;
                        v3 = sub_hash(v3, b[index]);
                        index++;
                        v4 = sub_hash(v4, b[index]);
                        index++;
                    } while (index <= limit);

                    h32 = rotate_left(v1, 1) + rotate_left(v2, 7) + rotate_left(v3, 12) + rotate_left(v4, 18);
                }
                else
                {
                    h32 = SEED + PRIME32_5;
                }

                h32 += (uint32_t)len * 4;

                while (index < len)
                {
                    h32 += b[index] * PRIME32_3;
                    h32 = rotate_left(h32, 17) * PRIME32_4;
                    index++;
                }

                h32 ^= h32 >> 15;
                h32 *= PRIME32_2;
                h32 ^= h32 >> 13;
                h32 *= PRIME32_3;
                h32 ^= h32 >> 16;

                return h32;
            }
        };

        struct Custom
        {
            inline static HashValue hash(void *buf)
            {
                const uint32_t _p1 = 15953071;
                const uint32_t _p2 = 37953119;
                const uint32_t _p3 = 73856093;
                const uint32_t _p4 = 93856897;
                uint32_t *b = (uint32_t *)buf;
                return ((_p1 * *(b + 0)) ^ (_p2 * *(b + 1)) ^ (_p3 * *(b + 2)) ^ (_p4 * *(b + 3)));
            }
        };

        struct Knuth
        {
            inline static HashValue hash(void *buf)
            {
                uint64_t *b = (uint64_t *)buf;
                // bitshift (middle bits contain more entropy) 8 instead of 16 to preserve lower bits for later range reduction
                return ((b[0] ^ b[1]) * 2654435761 >> 8);
            }
        };

    } // namespace hashing

    namespace reduction
    {

        struct Mod
        {
            inline static HashValue reduce(HashValue hash, uint32_t buckets)
            {
                return hash % buckets;
            }
        };

        struct FastRange
        {
            inline static HashValue reduce(HashValue hash, uint32_t buckets)
            {
                return ((uint64_t)hash * (uint64_t)buckets) >> 32;
            }
        };

        struct Identity
        {
            inline static HashValue reduce(HashValue hash, uint32_t buckets)
            {
                return hash;
            }
        };

    } // namespace reduction

    using real = double;

    template <typename Value, typename HashFunction = hashing::Murmur, typename ReduceFunction = reduction::FastRange, size_t REHASH_ROUNDS = 5>
    class SpatialHash
    {
    public:
        struct HashBucket
        {
            int x;
            int y;
            unsigned int last_claimed = -1;
            std::vector<Value> data;
        };

        SpatialHash();

        SpatialHash(
            real cell_size,
            unsigned int table_size);

        void reset(real cell_size, unsigned int table_size);

        void insert_at_cell(
            int x,
            int y,
            Value &value,
            int salt = 0);

        void insert_at_point(
            real x,
            real y,
            Value &value,
            int salt = 0);

        void insert_at_aabb(
            real top_left_x,
            real top_left_y,
            real bottom_right_x,
            real bottom_right_y,
            Value &value,
            int salt = 0);

        void insert_at_segment(
            real start_x_coord,
            real start_y_coord,
            real end_x_coord,
            real end_y_coord,
            Value &value,
            int salt = 0);

        void query_at_cell(
            std::vector<Value> &result,
            int x,
            int y,
            int salt = 0);

        void query_at_point(
            std::vector<Value> &result,
            real x,
            real y,
            int salt = 0);

        void query_at_aabb(
            std::vector<Value> &result,
            real top_left_x,
            real top_left_y,
            real bottom_right_x,
            real bottom_right_y,
            int salt = 0);

        void query_at_segment(
            std::vector<Value> &result,
            real start_x_coord,
            real start_y_coord,
            real end_x_coord,
            real end_y_coord,
            int salt = 0);

        //private:
        int cell(real coordinate) const;
        HashBucket *get_bucket(int x, int y, int salt);

    private:
        real _inv_cell_size;
        unsigned int _table_size;
        unsigned int _current_round = 0; // pepper

        std::vector<HashBucket> _hash_table;
    };

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::SpatialHash()
        : _inv_cell_size(1.0 / 1.0),
          _table_size(1024)
    {
        _hash_table.resize(_table_size, HashBucket());
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::SpatialHash(
        real cell_size,
        unsigned int table_size)
        : _inv_cell_size(1.0 / cell_size),
          _table_size(table_size)
    {
        _hash_table.resize(_table_size, HashBucket());
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::reset(real cell_size, unsigned int table_size)
    {
        _inv_cell_size = 1.0 / cell_size;
        _current_round += 1;

        if (_table_size != table_size)
        {
            _table_size = table_size;
            _hash_table.resize(_table_size, HashBucket());
        }
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::insert_at_cell(
        int x,
        int y,
        Value &value,
        int salt)
    {
        HashBucket *bucket = get_bucket(x, y, salt);
        bucket->data.push_back(value);
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::insert_at_point(
        real x,
        real y,
        Value &value,
        int salt)
    {
        int cell_x = cell(x);
        int cell_y = cell(y);

        insert_at_cell(cell_x, cell_y, value, salt);
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::insert_at_aabb(
        real top_left_x,
        real top_left_y,
        real bottom_right_x,
        real bottom_right_y,
        Value &value,
        int salt)
    {
        int top_left_x_cell = cell(top_left_x);
        int top_left_y_cell = cell(top_left_y);
        int bottom_right_x_cell = cell(bottom_right_x);
        int bottom_right_y_cell = cell(bottom_right_y);

        for (int i = top_left_x_cell; i <= bottom_right_x_cell; i++)
        {
            for (int j = top_left_y_cell; j <= bottom_right_y_cell; j++)
            {
                insert_at_cell(i, j, value, salt);
            }
        }
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::insert_at_segment(
        real start_x_coord,
        real start_y_coord,
        real end_x_coord,
        real end_y_coord,
        Value &value,
        int salt)
    {
        int x0 = cell(start_x_coord);
        int y0 = cell(start_y_coord);
        int x1 = cell(end_x_coord);
        int y1 = cell(end_y_coord);

        int dx = std::abs(x1 - x0);
        int dy = -std::abs(y1 - y0);
        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;
        int err = dx + dy;
        int e2;

        while (true)
        {
            insert_at_cell(x0, y0, value, salt);

            if (x0 == x1 && y0 == y1)
                break;

            e2 = 2 * err;

            if (e2 > dy)
            {
                err += dy;
                x0 += sx;
            }
            else if (e2 < dx)
            {
                err += dx;
                y0 += sy;
            }
        }
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::query_at_cell(
        std::vector<Value> &result,
        int x,
        int y,
        int salt)
    {
        HashBucket *bucket = get_bucket(x, y, salt);
        result.insert(result.begin(), bucket->data.begin(), bucket->data.end());
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::query_at_point(
        std::vector<Value> &result,
        real x,
        real y,
        int salt)
    {
        int cell_x = cell(x);
        int cell_y = cell(y);

        query_at_cell(result, cell_x, cell_y, salt);
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::query_at_aabb(
        std::vector<Value> &result,
        real top_left_x,
        real top_left_y,
        real bottom_right_x,
        real bottom_right_y,
        int salt)
    {
        int top_left_x_cell = cell(top_left_x);
        int top_left_y_cell = cell(top_left_y);
        int bottom_right_x_cell = cell(bottom_right_x);
        int bottom_right_y_cell = cell(bottom_right_y);

        for (int i = top_left_x_cell; i <= bottom_right_x_cell; i++)
        {
            for (int j = top_left_y_cell; j <= bottom_right_y_cell; j++)
            {
                query_at_cell(result, i, j, salt);
            }
        }
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    void SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::query_at_segment(
        std::vector<Value> &result,
        real start_x_coord,
        real start_y_coord,
        real end_x_coord,
        real end_y_coord,
        int salt)
    {
        int x0 = cell(start_x_coord);
        int y0 = cell(start_y_coord);
        int x1 = cell(end_x_coord);
        int y1 = cell(end_y_coord);

        int dx = std::abs(x1 - x0);
        int dy = std::abs(y1 - y0);
        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;
        int err = dx + dy;
        int e2;

        while (true)
        {
            query_at_cell(result, x0, y0, salt);

            if (x0 == x1 && y0 == y1)
                break;

            e2 = 2 * err;

            if (e2 > dy)
            {
                err += dy;
                x0 += sx;
            }
            else if (e2 < dx)
            {
                err += dx;
                y0 += sy;
            }
        }
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    int SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::cell(real coordinate) const
    {
        return std::floor(coordinate * _inv_cell_size);
    }

    template <typename Value, typename HashFunction, typename ReduceFunction, size_t REHASH_ROUNDS>
    typename SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::HashBucket *
    SpatialHash<Value, HashFunction, ReduceFunction, REHASH_ROUNDS>::get_bucket(int x, int y, int salt)
    {
        uint32_t buf[4];
        buf[0] = x;
        buf[1] = y;
        buf[2] = salt;
        buf[3] = _current_round;

        HashBucket *bucket;

        for (size_t i = 0; i <= REHASH_ROUNDS; i++)
        {
            buf[3] += 1;
            bucket = &_hash_table[ReduceFunction::reduce(HashFunction::hash((void *)buf), _table_size)];

            if (bucket->last_claimed != _current_round)
            {
                bucket->last_claimed = _current_round;
                bucket->data.clear();
                bucket->x = x;
                bucket->y = y;
                break;
            }
            else if (bucket->x == x && bucket->y == y)
            {
                break;
            }
        }

        return bucket;
    }

}; // namespace shash

#endif
