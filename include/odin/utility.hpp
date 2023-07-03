// Define macro for changing tbb::concurrent_vector to std::vector
#define CONVERT_TO_STD_VECTOR(name) std::vector<decltype(name)::value_type>(name.begin(), name.end())