import gzip
import logging
import pathlib
import shutil
from dataclasses import dataclass, field
from ftplib import FTP  # nosec
from io import StringIO
from typing import List, Set

from shared.utils.crypto import CryptoUtils
from shared.utils.path import PathUtils
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.enums import SyncItemStatus
from substances_core.domain.syncs.services.extractor import SyncExtractor

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PubchemExtractorOptions:
    """PubchemFetcherOptions"""

    host: str
    """URI to the provider data source."""

    data_location: str
    """Path to the sources files."""

    data_dir: str
    """Path to the directory used to store extracted data."""

    items_limit: int = 1
    """Maximum items to retrieve if no index specified."""

    indexes: Set[str] = field(default_factory=set)
    """Define specific files to extract."""


class PuchemSubstanceExtractor(SyncExtractor):
    """
    PuchemSubstanceExtractor
    """

    def __init__(self, options: PubchemExtractorOptions) -> None:
        self._options = options
        self._ftp_client = FTP(options.host)  # nosec

    def get_data_location(self) -> str:
        return self._options.data_location

    def extract(self) -> List[SyncItem]:
        """
        Performs extraction of SDF archives stored on PubChem
        FTP server.
        """
        self._ftp_client.login()

        # Try navigating to datasource locaion
        try:
            self._ftp_client.cwd(self.get_data_location())
        except Exception as err:
            LOGGER.error(err)
            return []

        items = self.__get_extractable_items(self._options.indexes)
        for item in items:
            try:
                self.__download_archive(item)
                self.__extract_archive(item)
                item.status = SyncItemStatus.EXTRACTED
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.error(err)
                item.status = SyncItemStatus.ERROR
        return items

    # Step 1
    def __get_extractable_items(self, indexes: Set[str] = set()) -> List[SyncItem]:
        """Extract list of archives containing SDF file.

        Returns:
            List[SyncItem]: list of archives available for synchronization.
        """
        files: List[str] = []
        # pylint: disable=unnecessary-lambda
        self._ftp_client.retrlines("NLST", lambda file: files.append(file))

        if len(indexes) > 0:
            files = [file for file in files if file in indexes]

        return self.__generate_items(
            [file for file in files if ".md5" not in file and "README" not in file]
        )

    # Step 1.1
    def __generate_items(self, files: List[str]) -> List[SyncItem]:
        items: List[SyncItem] = []

        max_items = (
            len(files)
            if self._options.items_limit > len(files)
            else self._options.items_limit
        )

        for file in files[0:max_items]:
            try:
                chain: str = self.__extract_checksum(f"{file}.md5")
                checksum = chain.split("  ")[0]
                # FTP use posix path
                location = pathlib.PurePosixPath.joinpath(
                    pathlib.PurePosixPath(self.get_data_location()), file
                )
                items.append(
                    SyncItem(token=checksum, location=str(location), name=file)
                )
                LOGGER.info("Checkin pubchem SDF archive: %r", chain)
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.error(err)
        return items

    def __extract_checksum(self, location: str) -> str:
        checkum_file = pathlib.PurePosixPath(location)
        try:
            reader = StringIO()
            self._ftp_client.retrlines("RETR " + checkum_file.name, reader.write)
            return reader.getvalue()
        except Exception as err:
            LOGGER.error(err)
            raise

    # Step 2
    def __download_archive(self, item: SyncItem) -> None:
        # archive = pathlib.PurePosixPath(item.location)
        if item.name is None:
            raise ValueError

        try:
            output_path = PathUtils.join_path(self._options.data_dir, item.name)
            # self._ftp_client.cwd(str(archive.parents[0]))
            self._ftp_client.retrbinary(
                "RETR " + item.name, open(str(output_path), "wb").write
            )
        except Exception:  # pylint: disable=broad-except
            raise

    # Step 3
    def __extract_archive(self, item: SyncItem) -> None:
        if item.name is None:
            raise ValueError

        archive_input_path = PathUtils.join_path(self._options.data_dir, item.name)
        archive_output_path = PathUtils.join_path(self._options.data_dir, item.token)
        try:
            if self.__validate_checksum(archive_input_path, item.token):
                with gzip.open(archive_input_path, "r") as gzip_input:
                    with open(archive_output_path, "wb") as f_out:
                        shutil.copyfileobj(gzip_input, f_out)
            else:
                raise RuntimeError
        except Exception:  # pylint: disable=broad-except
            raise

    def __validate_checksum(self, archive_path: str, checksum: str) -> bool:
        try:
            if CryptoUtils.md5sum_file(archive_path).casefold() == checksum.casefold():
                return True
            else:
                return False
        except Exception:  # pylint: disable=broad-except
            return False
