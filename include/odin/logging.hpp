
#ifndef LOGGING_HPP
#define LOGGING_HPP

#define ANSI_COLOR_RESET "\x1b[0m"
#define ANSI_COLOR_BLACK "\x1b[30m"
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_WHITE "\x1b[37m"

// Utility functions to colorize your log messages
std::string Colorize(const std::string &message, const std::string &color_code) {
    return color_code + message + ANSI_COLOR_RESET;
}

std::string Red(const std::string &message) {
    return Colorize(message, ANSI_COLOR_RED);
}

std::string Green(const std::string &message) {
    return Colorize(message, ANSI_COLOR_GREEN);
}

std::string Blue(const std::string &message) {
    return Colorize(message, ANSI_COLOR_BLUE);
}


#ifdef ODIN_USE_GLOG
#include <glog/logging.h>
#include <iomanip>
using namespace google;

//#ifdef GLOG_CUSTOM_PREFIX_SUPPORT
//struct LogMessageInfo {
//    explicit LogMessageInfo(const char* const severity_,
//                            const char* const filename_,
//                            const int& line_number_,
//                            const int& thread_id_,
//                            const LogMessageTime& time_):
//                                                           severity(severity_), filename(filename_), line_number(line_number_),
//                                                           thread_id(thread_id_), time(time_)
//    {}
//
//    const char* const severity;
//    const char* const filename;
//    const int &line_number;
//    const int &thread_id;
//    const LogMessageTime& time;
//};

enum Justification { LEFT,
                     RIGHT,
                     CENTER };

std::string PadString(const std::string &input, size_t length, char pad_char, Justification justification) {
    if (input.size() >= length) {
        return input;
    }

    std::string padding(length - input.size(), pad_char);

    switch (justification) {
        case Justification::LEFT:
            return padding + input;
        case Justification::RIGHT:
            return input + padding;
        case Justification::CENTER: {
            size_t pad_left  = padding.size() / 2;
            size_t pad_right = padding.size() - pad_left;
            return padding.substr(0, pad_left) + input + padding.substr(0, pad_right);
        }
        default:
            return input;// Return original string as default case.
    }
}

#ifdef GLOG_CUSTOM_PREFIX_SUPPORT
void CustomPrefix(std::ostream &stream, const google::LogMessageInfo &info, void *) {
    const char *severity_color = ANSI_COLOR_RESET;
    std::string severity_padded;

    if (strcmp(info.severity, "INFO") == 0) {
        severity_color  = ANSI_COLOR_BLUE;
        severity_padded = PadString("INFO", 5, ' ', Justification::RIGHT);
    } else if (strcmp(info.severity, "WARNING") == 0) {
        severity_color  = ANSI_COLOR_YELLOW;
        severity_padded = PadString("WARN", 5, ' ', Justification::RIGHT);
    } else if (strcmp(info.severity, "ERROR") == 0) {
        severity_color  = ANSI_COLOR_RED;
        severity_padded = PadString("ERR!", 5, ' ', Justification::RIGHT);
    } else if (strcmp(info.severity, "VERBOSE") == 0) {
        severity_color  = ANSI_COLOR_MAGENTA;
        severity_padded = PadString("VLOG", 5, ' ', Justification::RIGHT);
    } else if (strcmp(info.severity, "FATAL") == 0) {
        severity_color  = ANSI_COLOR_CYAN;
        severity_padded = PadString("FATAL", 5, ' ', Justification::RIGHT);
    } else {
        severity_color  = ANSI_COLOR_WHITE;
        severity_padded = PadString("UNKNOWN", 5, ' ', Justification::RIGHT);
    }

    stream << '[' << " " << severity_color << severity_padded << ANSI_COLOR_RESET << " | "
           << std::setw(4) << 1900 + info.time.year() << '-'
           << std::setw(2) << 1 + info.time.month() << '-'
           << std::setw(2) << info.time.day() << 'T'
           << std::setw(2) << info.time.hour() << ':'
           << std::setw(2) << info.time.min() << ':'
           << std::setw(2) << info.time.sec()
           //           << "."
           //           <<
           //            std::setw(6) << info.time.usec()
           << " | "
           << std::setfill(' ') << std::setw(5)
           << info.thread_id << std::setfill('0') << " | "
           << info.filename << ':' << info.line_number << " ] ";
}
#endif// GLOG_CUSTOM_PREFIX_SUPPORT

#define ODIN_IF_LOGGING_ENABLED(x) x
#define ODIN_IF_VERBOSITY_MATCHES(Expr, Level) \
    ((Level) <= FLAGS_v ? (Expr) : "")
#define ODIN_LOG_INFO LOG(INFO)
#define ODIN_LOG_WARNING LOG(WARNING)
#define ODIN_LOG_ERROR LOG(ERROR)
#define ODIN_LOG_FATAL LOG(FATAL)
#define ODIN_VLOG(verbosity) VLOG(verbosity)
#define ODIN_SET_VLOG_LEVEL(LEVEL) (FLAGS_v = LEVEL)

#ifdef GLOG_CUSTOM_PREFIX_SUPPORT
#define GOOGLE_INIT_LOGGING(NAME) \
    google::InitGoogleLogging(NAME, CustomPrefix)
#else
#define GOOGLE_INIT_LOGGING(NAME) \
    google::InitGoogleLogging(NAME)
#endif

#define ODIN_SET_LOG_DESTINATION(LEVEL, DESTINATION) \
    google::SetLogDestination(LEVEL, DESTINATION);

//#define INFO google::GLOG_INFO
//#define WARNING google::GLOG_WARNING
//#define ERROR google::GLOG_ERROR
//#define FATAL google::GLOG_FATAL


#define INIT_ODIN_LOGGING(NAME, LOG_DIR)                       \
    do {                                                       \
        GOOGLE_INIT_LOGGING(NAME);                             \
        google::SetLogDestination(google::GLOG_INFO, LOG_DIR); \
        FLAGS_logbufsecs                = 0;                   \
        FLAGS_colorlogtostderr          = true;                \
        FLAGS_log_prefix                = true;                \
        FLAGS_alsologtostderr           = true;                \
        FLAGS_max_log_size              = 100;                 \
        FLAGS_stop_logging_if_full_disk = true;                \
        FLAGS_logtostderr               = true;                \
        FLAGS_stderrthreshold           = 0;                   \
    } while (0)

#else
class NullStream {
public:
    template<typename T>
    NullStream &operator<<(const T &) { return *this; }
};

static NullStream nullstream;
#define ODIN_LOG_INFO nullstream
#define ODIN_LOG_WARNING nullstream
#define ODIN_LOG_ERROR nullstream
#define ODIN_LOG_FATAL nullstream
#define ODIN_VLOG(verbosity) nullstream
#define ODIN_SET_VLOG_LEVEL(LEVEL) (void) 0
#define INIT_ODIN_LOGGING(NAME, LOG_DIR) (void) 0
#define ODIN_IF_LOGGING_ENABLED(x) ((void) 0)
#define ODIN_IF_VERBOSITY_MATCHES(verbosity, x) ""
#define ODIN_SET_LOG_DESTINATION(LEVEL, DESTINATION) (void) 0

#endif// ODIN_USE_GLOG


#endif// LOGGING_HPP