import logging
import os
import pathlib
import sys

LOGGER = logging.getLogger("AppRoot")

sys.path.insert(
    1, str(pathlib.Path.joinpath(pathlib.Path(os.path.dirname(__file__)), "grpc"))
)

LOGGER.info("sys.path: %r", "; ".join(sys.path))
