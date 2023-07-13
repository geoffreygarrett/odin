
#ifndef RANDOM_HPP
#define RANDOM_HPP
//https://stackoverflow.com/questions/65871948/same-random-numbers-in-c-as-computed-by-python3-numpy-random-rand
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

namespace odin {

    class np_mt19937 {
    private:
        std::mt19937 engine;

    public:
        explicit np_mt19937(unsigned seed = std::random_device{}()) : engine(seed) {}

        double rand() {
            std::uint32_t a = engine() >> 5;
            std::uint32_t b = engine() >> 6;
            return (a * 67108864.0 + b) / 9007199254740992.0;
        }

        unsigned int randint() {
            return engine();
        }

        void save(const std::string &filepath) const {
            std::ofstream               os(filepath, std::ios::binary);
            cereal::BinaryOutputArchive archive(os);
            archive(*this);
        }

        void load(const std::string &filepath) {
            std::ifstream              is(filepath, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(*this);
        }

        // Serialization function definition
        template<class Archive>
        void serialize(Archive &ar) {
            ar(engine);
        }
    };

}// namespace odin

#endif//RANDOM_HPP
