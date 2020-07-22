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

#include <iostream>
#include <algorithm>
#include <vector>
#include <chrono>
#include "SpatialHash.h"

class SpatialHashTest
{
    using Id = uint32_t;
    using real = float;

    struct Object
    {
        real x;
        real y;
        int category;
        Id value;
    };

public:
    SpatialHashTest(int test_size, std::vector<real> load_factors, real cell_size, real world_size)
        : _test_size(test_size),
          _load_factors(load_factors),
          _cell_size(cell_size),
          _world_size(world_size)
    {
        init_test_data();
    }

    int test_all()
    {
        int result = 1;

        result *= test_get_cell();
        result *= test_insert_query_point();
        result *= test_insert_query_aabb();
        result *= test_insert_query_segment();

        return result;
    }

    int test_insert_query_segment()
    {
        std::cout << "Test insert query segment" << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        shash::SpatialHash<Id, shash::hashing::Murmur, shash::reduction::FastRange, 10> spatial_hash(1.0, 1000);

        Id val1 = 1;
        Id val2 = 2;

        spatial_hash.insert_at_segment(0.0, 0.0, 20.0, 20.0, val1, 1);
        spatial_hash.insert_at_segment(10.0, 0.0, 0.0, 30.0, val2, 1);

        // TODO

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "\tDuration: " << duration << " milliseconds " << std::endl;

        return 1;
    }

    int test_insert_query_aabb()
    {
        std::cout << "Test insert query aabb" << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();

        shash::SpatialHash<Id, shash::hashing::Murmur, shash::reduction::FastRange, 10> spatial_hash(1.0, 1000);

        Id val1 = 1;
        Id val2 = 2;

        spatial_hash.insert_at_aabb(0.0, 0.0, 20.0, 20.0, val1, 1);
        spatial_hash.insert_at_aabb(10.0, 10.0, 30.0, 30.0, val2, 1);

        std::vector<Id> result;
        spatial_hash.query_at_aabb(result, 18.0, 18.0, 20.0, 20.0, 1);
        if (result.size() != 18)
        {
            std::cout << "\tFAIL, did not find some data!!" << std::endl;
            return 0;
        }

        result.clear();
        spatial_hash.query_at_point(result, 20.0, 20.0, 1);
        if (result.size() != 2)
        {
            std::cout << "\tFAIL, did not find some data!!" << std::endl;
            return 0;
        }

        result.clear();
        spatial_hash.query_at_point(result, 1.0, 1.0, 1);
        if (result.size() != 1)
        {
            std::cout << "\tFAIL, did not find some data!!" << std::endl;
            return 0;
        }

        result.clear();
        spatial_hash.query_at_point(result, 25.0, 25.0, 1);
        if (result.size() != 1)
        {
            std::cout << "\tFAIL, did not find some data!!" << std::endl;
            return 0;
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "\tDuration: " << duration << " milliseconds " << std::endl;

        return 1;
    }

    int test_insert_query_point()
    {
        std::cout << "Test insert query point" << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();

        shash::SpatialHash<Id, shash::hashing::Murmur, shash::reduction::FastRange, 5> spatial_hash(_cell_size, 42);

        for (auto &load : _load_factors)
        {
            spatial_hash.reset(_cell_size, (int)(_test_size / load));
            real collisions = 0.0;

            for (auto &e : _test_data)
            {
                spatial_hash.insert_at_point(e.x, e.y, e.value, e.category);
            }

            std::vector<Id> result;
            for (auto &e : _test_data)
            {
                result.clear();
                spatial_hash.query_at_point(result, e.x, e.y, e.category);

                collisions += (result.size() - std::ceil(load)) / result.size();

                if (std::find(result.begin(), result.end(), e.value) == result.end())
                {
                    std::cout << "\tFAIL, did not find some data!!" << std::endl;
                    return 0;
                }
            }

            std::cout
                << "\tLoad: " << load
                << " \tCollisions: " << (real)collisions / _test_size
                << std::endl;
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "\tDuration: " << duration << " milliseconds " << std::endl;

        return 1;
    }

    int test_get_cell()
    {
        std::cout << "Test get_cell / hash collisions" << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        shash::SpatialHash<Id, shash::hashing::Murmur, shash::reduction::FastRange, 5> spatial_hash(_cell_size, 42);

        for (auto &load : _load_factors)
        {
            spatial_hash.reset(_cell_size, (int)(_test_size / load));
            int collisions = 0;

            for (auto &e : _test_data)
            {
                int cell_x = spatial_hash.cell(e.x);
                int cell_y = spatial_hash.cell(e.y);
                auto bucket = spatial_hash.get_bucket(cell_x, cell_y, e.category);
                bucket->data.push_back(e.value);
            }

            for (auto &e : _test_data)
            {
                int cell_x = spatial_hash.cell(e.x);
                int cell_y = spatial_hash.cell(e.y);
                auto bucket = spatial_hash.get_bucket(cell_x, cell_y, e.category);

                if (bucket->x != cell_x || bucket->y != cell_y)
                    collisions++;

                if (std::find(bucket->data.begin(), bucket->data.end(), e.value) == bucket->data.end())
                {
                    std::cout << "\tFAIL, did not find some data!!" << std::endl;
                    return 0;
                }
            }

            std::cout
                << "\tLoad: " << load
                << " \tCollisions: " << (real)collisions / _test_size
                << std::endl;
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "\tDuration: " << duration << " milliseconds " << std::endl;

        return 1;
    }

private:
    void init_test_data()
    {
        for (size_t i = 0; i < _test_size; i++)
        {
            Id entry = (Id)std::rand();
            real x = (-_world_size) + static_cast<real>(rand()) / (static_cast<real>(RAND_MAX / (_world_size * 2)));
            real y = (-_world_size) + static_cast<real>(rand()) / (static_cast<real>(RAND_MAX / (_world_size * 2)));
            int category = rand() % 256;

            _test_data.push_back({x, y, category, entry});
        }
    }

private:
    int _test_size;
    std::vector<real> _load_factors;
    real _cell_size;
    real _world_size;
    std::vector<Object> _test_data;
};

int main()
{
    SpatialHashTest test(100000, {0.1, 0.3, 0.5, 0.7, 1.0, 2.0}, 1.0, 1000000.0);
    if (!test.test_all())
        std::cout << "All Tests FAILED" << std::endl;
    else
        std::cout << "All Tests SUCCEEDED" << std::endl;

    return 0;
}