import hashlib


class CryptoUtils:
    """Utils class with cryptographic helpers method."""
    @staticmethod
    def md5sum_file(filepath: str, blocksize=65536):
        md5 = hashlib.md5(usedforsecurity=False)
        with open(filepath, "rb") as file:
            for block in iter(lambda: file.read(blocksize), b""):
                md5.update(block)
        return md5.hexdigest()
