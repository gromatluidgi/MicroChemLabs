import logging
from enum import Enum

from .config import settings

LOG_FORMAT_DEBUG = "%(levelname)s:%(message)s:%(pathname)s:%(funcName)s:%(lineno)d"


class LogLevels(str, Enum):
    info = "INFO"
    warn = "WARN"
    error = "ERROR"
    debug = "DEBUG"


class LoggingFilter(logging.Filter):
    """LoggingFilter"""

    def filter(self, record):
        """Determine which log records to output.
        Returns 0 for no, nonzero for yes.
        """
        if record.funcName.startswith("pika."):
            return False
        return True


def configure_logging():
    logging.getLogger().addFilter(LoggingFilter())

    log_level = str(settings.LOG_LEVEL).upper()  # cast to string
    log_levels = [level for level in LogLevels]

    if log_level not in log_levels:
        # we use error as the default log level
        logging.basicConfig(level=LogLevels.error)
        return

    if log_level == LogLevels.debug:
        logging.basicConfig(level=log_level, format=LOG_FORMAT_DEBUG)
        return

    logging.basicConfig(level=log_level)
