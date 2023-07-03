//#ifndef PROTOBUF_MEMORY_POLICY_H
//#define PROTOBUF_MEMORY_POLICY_H
//
//#include <string>
//#include <vector>
//#include "memory_policy.hpp"
//#include <core/policy/serialization/serialization_policy.hpp>
//
//// Delimiter character used to separate serialized data items in a string.
//// Make sure this character does not appear in the serialized data.
//const char kDelimiterCharacter = '\n';
//
//template<typename Data, typename Container, typename Serializer = SerializationPolicy<Data>>
//class ProtobufMemoryPolicy : public MemoryPolicy<ProtobufMemoryPolicy<Data, Container, Serializer>, Container> {
//public:
//    void store_impl(const Container &data) {
//        // Serialize each data item in the container into a string.
//        Serializer serializer;
//        std::string serialized_data;
//        for (const auto &item: data) {
//            serialized_data += serializer.serialize(item);
//            serialized_data += kDelimiterCharacter;
//        }
//
//        // TODO: Store the serialized data somewhere...
//    }
//
//    Container retrieve_impl() {
//        // TODO: Retrieve the serialized data from somewhere...
//        std::string serialized_data;
//
//        // Parse each item from the serialized data string and convert back to a container.
//        Container data;
//        while (!serialized_data.empty()) {
//            size_t delimiter_pos = serialized_data.find(kDelimiterCharacter);
//            std::string single_item_serialized_data = serialized_data.substr(0, delimiter_pos);
//            serialized_data = serialized_data.substr(delimiter_pos + 1);
//
//            Data item = serializer.deserialize(single_item_serialized_data);
//            data.push_back(item);
//        }
//        return data;
//    }
//
//    void delete_data_impl() {
//        // TODO: Delete the stored serialized data...
//    }
//};
//
//#endif // PROTOBUF_MEMORY_POLICY_H
