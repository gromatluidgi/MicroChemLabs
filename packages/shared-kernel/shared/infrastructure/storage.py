from abc import ABC, abstractmethod
from typing import IO


class Storage(ABC):
    @abstractmethod
    def copy_file_from_stream(self, stream: IO, target_path: str):
        """Copies data from a stream to a file in the storage.

        Args:
            stream (IO): https://docs.python.org/3/library/typing.html#typing.IO
            target_path (str):

        Raises:
            NotImplementedError: _description_
        """
        raise NotImplementedError

    def get_file_from_stream(self, stream):
        raise NotImplementedError


class FileStorage(Storage):
    ...
