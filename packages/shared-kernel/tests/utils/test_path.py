from shared.utils.path import PathUtils


def test_join_path():
    base_path = "toto/"
    others = "test.txt"

    result = PathUtils.join_path(base_path, others)

    assert result is not None
