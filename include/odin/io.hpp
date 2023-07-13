#ifndef ODIN_IO_HPP
#define ODIN_IO_HPP

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/common.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>
#include <sstream>
#include <string>

namespace odin {
    template<typename T>
    std::string to_json(const T &object) {
        std::ostringstream os;
        {
            cereal::JSONOutputArchive archive(os);
            archive(cereal::make_nvp("data", object));
        }
        return os.str();
    }

    template<typename T>
    void save_to_json(const T &object, const std::string &filename) {
        std::ofstream ofs(filename);
        if (ofs) {
            ofs << to_json(object);
            ofs.close();
        } else {
            throw std::runtime_error("Failed to open the file: " + filename);
        }
    }

    template<typename T>
    T from_json(const std::string &json_str) {
        std::istringstream is(json_str);
        T                  object;

        {
            cereal::JSONInputArchive archive(is);
            archive(cereal::make_nvp("data", object));
        }
        return object;
    }

    template<typename T>
    T load_from_json(const std::string &filename) {
        std::ifstream ifs(filename);
        if (!ifs) {
            throw std::runtime_error("Failed to open the file: " + filename);
        }
        return from_json<T>(std::string(std::istreambuf_iterator<char>(ifs),
                                        std::istreambuf_iterator<char>()));
    }

    template<typename T>
    std::string to_binary(const T &object, bool portable = false) {
        std::stringstream ss;
        if (portable) {
            cereal::PortableBinaryOutputArchive archive(ss);
            archive(cereal::make_nvp("data", object));
        } else {
            cereal::BinaryOutputArchive archive(ss);
            archive(cereal::make_nvp("data", object));
        }
        return ss.str();
    }

    template<typename T>
    void save_to_binary(const T &object, const std::string &filename, bool portable = false) {
        std::ofstream ofs(filename, std::ios::binary);
        if (ofs) {
            ofs << to_binary(object, portable);
            ofs.close();
        } else {
            throw std::runtime_error("Failed to open the file: " + filename);
        }
    }

    template<typename T>
    T from_binary(const std::string &binary_str, bool portable = false) {
        std::istringstream is(binary_str, std::ios::binary);
        T                  object;
        if (portable) {
            cereal::PortableBinaryInputArchive archive(is);
            archive(cereal::make_nvp("data", object));
        } else {
            cereal::BinaryInputArchive archive(is);
            archive(cereal::make_nvp("data", object));
        }
        return object;
    }

    template<typename T>
    T load_from_binary(const std::string &filename, bool portable = false) {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            throw std::runtime_error("Failed to open the file: " + filename);
        }
        std::string binary_str((std::istreambuf_iterator<char>(ifs)),
                               std::istreambuf_iterator<char>());
        return from_binary<T>(binary_str, portable);
    }


    template<typename T>
    std::string to_xml(const T &object) {
        std::ostringstream os;
        {
            cereal::XMLOutputArchive archive(os);
            archive(cereal::make_nvp("data", object));
        }
        return os.str();
    }

    template<typename T>
    void save_to_xml(const T &object, const std::string &filename) {
        std::ofstream ofs(filename);
        if (ofs) {
            ofs << to_xml(object);
            ofs.close();
        } else {
            throw std::runtime_error("Failed to open the file: " + filename);
        }
    }

    template<typename T>
    T from_xml(const std::string &xml_str) {
        std::istringstream is(xml_str);
        T                  object;
        {
            cereal::XMLInputArchive archive(is);
            archive(cereal::make_nvp("data", object));
        }
        return object;
    }

    template<typename T>
    T load_from_xml(const std::string &filename) {
        std::ifstream ifs(filename);
        if (!ifs) {
            throw std::runtime_error("Failed to open the file: " + filename);
        }
        return from_xml<T>(std::string(std::istreambuf_iterator<char>(ifs),
                                       std::istreambuf_iterator<char>()));
    }

}// namespace odin

#endif//ODIN_IO_HPP