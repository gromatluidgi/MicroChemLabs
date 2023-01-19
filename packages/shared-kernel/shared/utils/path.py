import pathlib
import tempfile


class PathUtils:
    """PathUtils"""

    @staticmethod
    def absolute_path(path: str) -> str:
        typed_path = pathlib.Path(path)
        if typed_path.is_file():
            return str(typed_path.parent.absolute())

        return str(typed_path.absolute())

    @staticmethod
    def join_path(base_path: str, *path: str) -> str:
        return str(pathlib.Path(base_path).joinpath(*path))

    @staticmethod
    def root_path() -> str:
        return PathUtils.absolute_path(pathlib.Path(pathlib.Path.home()).parts[0])

    @staticmethod
    def get_parent_dir(path: str, level=0) -> str:
        return str(pathlib.Path(path).parents[level])

    @staticmethod
    def generate_temp_dir(dirname: str) -> str:
        tmp_dir = pathlib.Path(tempfile.gettempdir()).joinpath(dirname)
        if not pathlib.Path.exists(tmp_dir):
            pathlib.Path.mkdir(tmp_dir)
        return str(tmp_dir)
