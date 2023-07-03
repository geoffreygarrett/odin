# Logging Policy

Effective logging is crucial to understand the behavior of our software and to diagnose problems. In
this project, we are using the Google Logging Library (glog) that provides a flexible and efficient
logging system. To ensure consistency and readability across the project, it's necessary to define a
logging policy. This policy document will define the guiding philosophy, requirements and severity
levels of the logging system.

## Guiding Philosophy

Our logging system aims to provide a precise understanding of the software's operational flow and
potential issues. We aim to make it as intuitive as possible to aid in debugging and system
comprehension. The guiding philosophy is "Informative, but not verbose" - each log should provide
clear, concise, and useful information.

## Requirements

The logging system should:

- **Informative**: Each log message should provide clear and concise information about what the
  system is doing.
- **Consistent**: Using the same terminology and style across all logs to minimize confusion.
- **Balanced**: The logging system should strike a balance between providing enough information for
  debugging purposes and avoiding information overload.

## Log Levels

We integrate Google's logging severity levels (INFO, WARNING, ERROR, and FATAL) with verbosity
levels in Google's VLOG. Here is our proposed hierarchy:

- **LOG(FATAL)**: For unrecoverable errors which require immediate program termination, e.g., "
  Cannot open crucial data file: termination required".
- **LOG(ERROR)**: For serious error events that might influence the program to abort, e.g., "
  Database connection lost".
- **VLOG(0)**: For important system-wide information and non-critical errors, e.g., "Service X
  started".
- **LOG(WARNING)**: For potentially harmful situations deserving users' attention, e.g., "Disk usage
  is high, might run out of space".
- **VLOG(10)**: For high-level operational messages indicating major system events, e.g., "User X
  has logged in".
- **LOG(INFO)**: For information representing the progress of the application, e.g., "Server has
  processed 1000 requests".
- **VLOG(20)** to **VLOG(30)**: For detailed operational information for specific subsystems,
  e.g., "Updating entity X in database".
- **VLOG(40)** to **VLOG(60)**: For detailed debugging information such as function entries/exits,
  variable values, e.g., "Entering loop with 100 iterations".
- **VLOG(70)** to **VLOG(90)**: For very detailed logs for diagnosing particularly difficult issues,
  e.g., "Output of complex data structure in raw format".

Please note, these log levels are not set in stone and should be adjusted as per the specific needs
and complexity of the system. The goal should always be to maximize the usefulness of the logs at
every level.

| Level               | Description                                                                    | Example                                               |
|---------------------|--------------------------------------------------------------------------------|-------------------------------------------------------|
| LOG(FATAL)          | Unrecoverable errors which require immediate program termination               | "Cannot open crucial data file: termination required" |
| LOG(ERROR)          | Serious error events that might influence the program to abort                 | "Database connection lost"                            |
| VLOG(0)             | Important system-wide information and non-critical errors                      | "Service X started"                                   |
| LOG(WARNING)        | Potentially harmful situations deserving users' attention                      | "Disk usage is high, might run out of space"          |
| VLOG(10)            | High-level operational messages indicating major system events                 | "User X has logged in"                                |
| LOG(INFO)           | Information representing the progress of the application                       | "Server has processed 1000 requests"                  |
| VLOG(20) - VLOG(30) | Detailed operational information for specific subsystems                       | "Updating entity X in database"                       |
| VLOG(40) - VLOG(60) | Detailed debugging information such as function entries/exits, variable values | "Entering loop with 100 iterations"                   |
| VLOG(70) - VLOG(90) | Very detailed logs for diagnosing particularly difficult issues                | "Output of complex data structure in raw format"      |

Remember that this table serves as a guideline and you can customize it to better suit your needs.
It is recommended to maintain consistency throughout your application for easier understanding and
troubleshooting.