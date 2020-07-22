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
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <array>
#include <cstdlib>
#include <chrono>
#include "SpatialHash.h"

struct Element
{
    int32_t x;
    int32_t y;
    uint32_t salt;
    uint32_t pepper;
};

int HASH_ROUNDS = 5;

static std::array<Element, 50000000> ELEMENTS;

void init_test_data()
{
    const int32_t min_coord = -100000;
    const int32_t max_coord = 100000;

    for (size_t i = 0; i < ELEMENTS.size(); i++)
    {
        ELEMENTS[i].x = min_coord + std::rand() % ((max_coord + 1) - min_coord);
        ELEMENTS[i].y = min_coord + std::rand() % ((max_coord + 1) - min_coord);
        ELEMENTS[i].salt = std::rand() % 256;
        ELEMENTS[i].pepper = 42;
    }
}

int realistic_speed_benchmark(shash::HashFuncPtr h, shash::ReduceFuncPtr r)
{
    uint32_t result = 0;
    const uint32_t buckets = 4096;
    auto hash_rounds = HASH_ROUNDS;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < ELEMENTS.size(); i++)
    {
        auto e = (&ELEMENTS[i]);
        auto e_ptr = (void *)e;

        for (size_t i = 0; i < hash_rounds; i++)
        {
            result += r(h(e_ptr), buckets);
            e->pepper += 1;
            e->salt += 1;
            e->x += 1;
            e->y += 1;
        }
        e->pepper -= hash_rounds;
        e->salt -= hash_rounds;
        e->x -= hash_rounds;
        e->y -= hash_rounds;
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << duration << " milliseconds " << std::endl;

    return result;
}

int realistic_speed_benchmark_reducers(shash::HashFuncPtr h)
{
    int result = 0;
    std::cout << "\t  No Reduce: \t\t";
    result += realistic_speed_benchmark(h, shash::reduction::Identity::reduce);
    std::cout << "\t  Mod Reduce: \t\t";
    result += realistic_speed_benchmark(h, shash::reduction::Mod::reduce);
    std::cout << "\t  Fast Range Reduce: \t";
    result += realistic_speed_benchmark(h, shash::reduction::FastRange::reduce);

    return result;
}

int realistic_speed_benchmark_all()
{
    int result = 0;

    std::cout << "Speed Benchmark: " << HASH_ROUNDS << " Rounds" << std::endl;

    std::cout << "\tCustom" << std::endl;
    result += realistic_speed_benchmark_reducers(shash::hashing::Custom::hash);

    std::cout << "\tCustom" << std::endl;
    result += realistic_speed_benchmark_reducers(shash::hashing::Custom::hash);

    std::cout << "\tKnuth" << std::endl;
    result += realistic_speed_benchmark_reducers(shash::hashing::Knuth::hash);

    std::cout << "\txxHash" << std::endl;
    result += realistic_speed_benchmark_reducers(shash::hashing::xxHash::hash);

    std::cout << "\tMurmur" << std::endl;
    result += realistic_speed_benchmark_reducers(shash::hashing::Murmur::hash);

    return result;
}

int realistic_quality_benchmark(shash::HashFuncPtr h, shash::ReduceFuncPtr r)
{
    int result = 0;

    std::vector<double> load_factors = {0.1, 0.3, 0.5, 0.7, 1, 2.0};
    std::vector<unsigned int> bucket_sizes = {512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};

    for (auto &load : load_factors)
    {
        std::vector<double> collision_rates;
        double min_collisions = 1.0;
        double max_collisions = 0.0;
        unsigned int min_cells;
        unsigned int max_cells;
        for (auto buckets : bucket_sizes)
        {
            std::unordered_map<unsigned int, int> occurences;
            for (int i = 0; i < buckets; i++)
            {
                occurences[i] = 0;
            }

            unsigned int cells = (unsigned int)buckets * load;
            unsigned int collisions = 0;

            for (size_t i = 0; i < cells; i++)
            {
                auto e = (&ELEMENTS[i]);
                auto e_ptr = (void *)e;
                auto hash_rounds = HASH_ROUNDS;

                auto hash = r(h(e_ptr), buckets);
                for (size_t i = 0; i < hash_rounds - 1; i++)
                {
                    if (occurences[hash] < std::ceil(load))
                        break;

                    e->pepper += (i + 1);
                    hash = r(h(e_ptr), buckets);
                    e->pepper -= (i + 1);
                }

                occurences[hash] += 1;

                if (occurences[hash] > 1)
                    collisions += 1;
            }

            double collision_rate = (double)collisions / cells;
            collision_rates.push_back(collision_rate);

            if (min_collisions > collision_rate)
            {
                min_collisions = collision_rate;
                min_cells = buckets;
            }
            if (max_collisions < collision_rate)
            {
                max_collisions = collision_rate;
                max_cells = buckets;
            }
        }

        double avg_collision_rate = 0.0;
        for (auto &cr : collision_rates)
        {
            avg_collision_rate += cr;
        }
        avg_collision_rate /= collision_rates.size();

        double variance = 0.0;
        for (auto &cr : collision_rates)
        {
            variance += (cr - avg_collision_rate) * (cr - avg_collision_rate);
        }

        variance /= collision_rates.size();
        double deviation = std::sqrt(variance);

        std::cout
            << "\t  Load: " << std::fixed << std::setw(11) << std::setprecision(6) << load
            << "\t  Avg.: " << std::fixed << std::setw(11) << std::setprecision(6) << avg_collision_rate
            << "\t  Deviation: " << std::fixed << std::setw(11) << std::setprecision(6) << deviation
            << "\t  Extremes: " << std::fixed << std::setw(11) << std::setprecision(6) << max_collisions - min_collisions
            << std::endl;
    }

    return result;
}

int realistic_quality_benchmark_reducers(shash::HashFuncPtr h)
{
    int result = 0;
    std::cout << "\t Mod Reduce \t\t" << std::endl;
    result += realistic_quality_benchmark(h, shash::reduction::Mod::reduce);
    std::cout << "\t Fast Range Reduce \t" << std::endl;
    result += realistic_quality_benchmark(h, shash::reduction::FastRange::reduce);

    return result;
}

int realistic_quality_benchmark_all()
{
    int result = 0;

    std::cout << "Quality Benchmark: " << HASH_ROUNDS << " Rounds" << std::endl;

    std::cout << "\tCustom" << std::endl;
    result += realistic_quality_benchmark_reducers(shash::hashing::Custom::hash);

    std::cout << "\tKnuth" << std::endl;
    result += realistic_quality_benchmark_reducers(shash::hashing::Knuth::hash);

    std::cout << "\txxHash" << std::endl;
    result += realistic_quality_benchmark_reducers(shash::hashing::xxHash::hash);

    std::cout << "\tMurmur" << std::endl;
    result += realistic_quality_benchmark_reducers(shash::hashing::Murmur::hash);

    return result;
}

int main()
{
    init_test_data();

    int result = 0;
    result += realistic_speed_benchmark_all();
    result += realistic_quality_benchmark_all();

    std::cout << result << std::endl;

    return 0;
}