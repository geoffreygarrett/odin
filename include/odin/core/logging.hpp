#ifdef ODIN_USE_GLOG
#include <glog/logging.h>
  #define ODIN_LOG_INFO(LOG) LOG(INFO) << LOG
  #define ODIN_LOG_WARNING(LOG) LOG(WARNING) << LOG
  #define ODIN_LOG_ERROR(LOG) LOG(ERROR) << LOG
#else
  #define ODIN_LOG_INFO(LOG) (void)0
  #define ODIN_LOG_WARNING(LOG) (void)0
  #define ODIN_LOG_ERROR(LOG) (void)0
#endif