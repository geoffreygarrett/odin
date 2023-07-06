#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <memory>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/async.h>

#define ANSI_COLOR_RESET "\x1b[0m"
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_WHITE "\x1b[37m"

enum LogLevel {
    TRACE = spdlog::level::trace,
    DEBUG = spdlog::level::debug,
    INFO = spdlog::level::info,
    WARN = spdlog::level::warn,
    ERROR = spdlog::level::err,
    CRITICAL = spdlog::level::critical
};

template<typename Derived>
class BaseLogger {
    std::once_flag async_init_flag;

public:
    std::shared_ptr<spdlog::logger> logger;

    BaseLogger() {
        std::call_once(async_init_flag, []() { spdlog::init_thread_pool(8192, 1); });
    }

    template<typename T>
    Derived &operator<<(const T &msg) {
        std::ostringstream stream;
        stream << msg;
        logger->info(stream.str());
        return static_cast<Derived &>(*this);
    }

    template<typename FormatString, typename... Args>
    void info(const FormatString &fmt, const Args &...args) {
        logger->info(fmt, args...);
    }

    template<typename FormatString, typename... Args>
    void debug(const FormatString &fmt, const Args &...args) {
        logger->debug(fmt, args...);
    }

    void set_level(spdlog::level::level_enum level) {
        logger->set_level(level);
    }

    template<typename FormatString, typename... Args>
    void warn(const FormatString &fmt, const Args &...args) {
        logger->warn(fmt, args...);
    }

    template<typename FormatString, typename... Args>
    void error(const FormatString &fmt, const Args &...args) {
        logger->error(fmt, args...);
    }

    template<typename FormatString, typename... Args>
    void critical(const FormatString &fmt, const Args &...args) {
        logger->critical(fmt, args...);
    }

    template<typename FormatString, typename... Args>
    void trace(const FormatString &fmt, const Args &...args) {
        logger->trace(fmt, args...);
    }
};

class ConsoleLogger : public BaseLogger<ConsoleLogger> {
public:
    ConsoleLogger() {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::info);
        logger = std::make_shared<spdlog::async_logger>("APP_LOGGER", console_sink, spdlog::thread_pool(),
                                                        spdlog::async_overflow_policy::block);
        logger->set_level(spdlog::level::trace);
        spdlog::register_logger(logger);
    }
};

class FileLogger : public BaseLogger<FileLogger> {
public:
    FileLogger(const std::string &filename) {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, true);
        file_sink->set_level(spdlog::level::trace);
        logger = std::make_shared<spdlog::async_logger>("APP_LOGGER", file_sink, spdlog::thread_pool(),
                                                        spdlog::async_overflow_policy::block);
        logger->set_level(spdlog::level::trace);
        spdlog::register_logger(logger);
    }
};

ConsoleLogger &global_console_logger() {
    static ConsoleLogger logger;
    return logger;
}

FileLogger &global_file_logger(const std::string &filename) {
    static FileLogger logger(filename);
    return logger;
}

class VLogger : public BaseLogger<VLogger> {
    int verbosity_level = 0;

public:
    VLogger(int v_level) {
        set_verbosity(v_level);
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::info);
        logger = std::make_shared<spdlog::async_logger>("APP_LOGGER", console_sink, spdlog::thread_pool(),
                                                        spdlog::async_overflow_policy::block);
        logger->set_level(spdlog::level::trace);
        spdlog::register_logger(logger);
    }

    void set_verbosity(int level) {
        verbosity_level = level;
    }

    int get_verbosity() const {
        return verbosity_level;
    }

    template<typename T>
    VLogger &operator<<(const T &msg) {
        if (get_verbosity() >= 20) {
            BaseLogger::operator<<(msg);
        }
        return *this;
    }
};

//#define ODIN_VLOG(v_level) VLogger(v_level) <<;

//#define ODIN_LOG_INFO StreamLogger(global_console_logger(), global_file_logger("logs/app.log"), spdlog::level::info)
//#define ODIN_LOG_WARNING StreamLogger(global_console_logger(), global_file_logger("logs/app.log"), spdlog::level::warn)
//#define ODIN_LOG_ERROR StreamLogger(global_console_logger(), global_file_logger("logs/app.log"), spdlog::level::err)
//#define ODIN_LOG_FATAL StreamLogger(global_console_logger(), global_file_logger("logs/app.log"), spdlog::level::critical)
//#define ODIN_LOG_TRACE StreamLogger(global_console_logger(), global_file_logger("logs/app.log"), spdlog::level::trace)

// Null stream that discards all input
class NullStream {
public:
    template<typename T>
    NullStream &operator<<(const T &) {
        return *this;
    }
};

#define ODIN_LOG_INFO NullStream()
#define ODIN_LOG_WARNING NullStream()
#define ODIN_LOG_ERROR NullStream()
#define ODIN_LOG_FATAL NullStream()
#define ODIN_LOG_TRACE NullStream()
#define ODIN_VLOG(v_level) NullStream()

//#define ODIN_LOG_INFO(...) global_console_logger().info(__VA_ARGS__); global_file_logger("logs/app.log").info(__VA_ARGS__)
//#define ODIN_LOG_WARNING(...) global_console_logger().warn(__VA_ARGS__); global_file_logger("logs/app.log").warn(__VA_ARGS__)
//#define ODIN_LOG_ERROR(...) global_console_logger().error(__VA_ARGS__); global_file_logger("logs/app.log").error(__VA_ARGS__)
//#define ODIN_LOG_FATAL(...) global_console_logger().critical(__VA_ARGS__); global_file_logger("logs/app.log").critical(__VA_ARGS__)
//#define ODIN_LOG_TRACE(...) global_console_logger().trace(__VA_ARGS__); global_file_logger("logs/app.log").trace(__VA_ARGS__)

#endif // LOGGING_HPP
